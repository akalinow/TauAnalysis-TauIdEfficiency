#include "TauAnalysis/TauIdEfficiency/plugins/ObjValEDNtupleProducer.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

typedef std::vector<unsigned> vunsigned;

namespace {
  // Function to build a nice name for the column
  std::string createBranchName(const std::string& ntupleName, const std::string& collection,
                               const std::string& varName, const int index=-1)
  {
    std::ostringstream out;
    out << ntupleName << "#" << collection;
    // Determine if this collection is indexed
    if ( index >= 0 ) {
      out << "@idx" << index;
    }
    out << "#" << varName;
    return out.str();
  }
}

void insertInto(std::vector<std::string>& a, const std::vector<std::string>& b)
{
  a.insert(a.begin(), b.begin(), b.end());
}

ObjValEDNtupleProducer::ObjValEDNtupleProducer(const edm::ParameterSet& cfg)
{
  std::cout << "<ObjValEDNtupleProducer::ObjValEDNtupleProducer>:" << std::endl;

  ntupleName_ = cfg.getParameter<std::string>("ntupleName");

  // Get list of ntuple sources to produce
  edm::ParameterSet cfgNtuples = cfg.getParameter<edm::ParameterSet>("sources");
  std::vector<std::string> ntupleNames = cfgNtuples.getParameterNamesForType<edm::ParameterSet>();

  for ( std::vector<std::string>::const_iterator ntupleName = ntupleNames.begin();
       ntupleName != ntupleNames.end(); ++ntupleName ) {
    std::cout << " reading configuration parameters for ntuple source = " << (*ntupleName) << std::endl;

    edm::ParameterSet cfgObjValExtractor = cfgNtuples.getParameter<edm::ParameterSet>(*ntupleName);
    std::string pluginTypeObjValExtractor = cfgObjValExtractor.getParameter<std::string>("pluginType");
    edm::InputTag cfgCollection = cfgObjValExtractor.getParameter<edm::InputTag>("src");

    bool isVectorExtractor = cfgObjValExtractor.getParameter<bool>("vector");
    
    // Determine if we want the whole collection or just specific indices. 
    // Default case: index = 0 only
    vunsigned indices;
    if ( cfgObjValExtractor.exists("indices") ) {
      indices = cfgObjValExtractor.getParameter<vunsigned>("indices");
    } else if ( isVectorExtractor ) {
      // if we want all objects, indicate it by leaving indices empty
    } else {
      // default case = take first index
      indices.push_back(0);
    }

    // Get the names of all the columns
    edm::ParameterSet columnsCfg = cfgObjValExtractor.getParameter<edm::ParameterSet>("columns");
    std::vector<std::string> columnNames;
    insertInto(columnNames, columnsCfg.getParameterNamesForType<std::string>());
    insertInto(columnNames, columnsCfg.getParameterNamesForType<edm::ParameterSet>());

    // Build each column source
    for( std::vector<std::string>::const_iterator columnName = columnNames.begin();
        columnName != columnNames.end(); ++columnName ) {
      edm::ParameterSet columnCfg;
      if ( columnsCfg.existsAs<std::string>(*columnName) ) {
	std::string columnValue = columnsCfg.getParameter<std::string>(*columnName);

        // concatentate relavant information into a single PSet to send to the ObjValExtractors
        columnCfg = cfgObjValExtractor;
	columnCfg.addParameter<std::string>("value", columnValue);
      } else if ( columnsCfg.existsAs<edm::ParameterSet>(*columnName) ) {
	columnCfg = columnsCfg.getParameter<edm::ParameterSet>(*columnName);
      } else assert(0);  

      // Check if we are taking specific indices
      if ( isVectorExtractor ) {
	// Build a vector getter
        ObjValVectorExtractorBase* objValExtractor =
            ObjValVectorExtractorPluginFactory::get()->create(pluginTypeObjValExtractor, columnCfg);
        std::string niceColumnName = createBranchName(ntupleName_, cfgCollection.label(), *columnName);
        ntupleVectorEntryType* ntupleEntry = new ntupleVectorEntryType(niceColumnName, objValExtractor, indices);
        ntupleVectorEntries_.push_back(ntupleEntry);
      } else {
        // Loop over each index
        for ( vunsigned::const_iterator index = indices.begin(); 
	      index != indices.end(); ++index ) {
          edm::ParameterSet columnCfgWithIndex = columnCfg;
          columnCfgWithIndex.addParameter<unsigned>("index", *index);
          // Build entry
          ObjValExtractorBase* objValExtractor =
              ObjValExtractorPluginFactory::get()->create(pluginTypeObjValExtractor, columnCfgWithIndex);
          std::string niceColumnName;
          if ( indices.size() == 1 && (*index) == 0 )
            niceColumnName = createBranchName(ntupleName_, cfgCollection.label(), *columnName);
          else
            niceColumnName = createBranchName(ntupleName_, cfgCollection.label(), *columnName, *index);
	  
          ntupleEntryType* ntupleEntry = new ntupleEntryType(niceColumnName, objValExtractor);
          ntupleEntries_.push_back(ntupleEntry);
        }
      } 
    }
  }

//--- add run, luminosity section and event numbers
  produces<edm::RunNumber_t>("run").setBranchAlias("run");
  produces<edm::LuminosityBlockNumber_t>("lumisection").setBranchAlias("lumisection");
  produces<edm::EventNumber_t>("event").setBranchAlias("event");

//--- register all the products with the framework
//    EK: count how many times we register each product
  std::map<std::string, int> name_counter;

  for ( std::vector<ntupleEntryType*>::const_iterator entry = ntupleEntries_.begin();
       entry != ntupleEntries_.end(); ++entry ) {
    std::string name = (*entry)->ntupleName_;
    name_counter[name] += 1;
    produces<double>(name).setBranchAlias(name);
  }

  for ( std::vector<ntupleVectorEntryType*>::const_iterator entry = ntupleVectorEntries_.begin();
       entry != ntupleVectorEntries_.end(); ++entry ) {
    std::string name = (*entry)->ntupleName_;
    name_counter[name] += 1;
    produces<std::vector<double> >(name).setBranchAlias(name);
  }

//--- make sure no variable is declared twice
  for ( std::map<std::string, int>::iterator imap = name_counter.begin();
	imap != name_counter.end(); ++imap ) {
    if ( imap->second > 1 ) {
      throw cms::Exception("DuplicateNtupleColumn") 
	<< "the ntuple variable: " << imap->first << " is declared " << imap->second << " times." 
	<< " It can only be declared once !! Please fix.";
    }
  }
  std::cout << "done." << std::endl;

  numEvents_processed_ = 0;
}

ObjValEDNtupleProducer::~ObjValEDNtupleProducer()
{
  std::cout << "<ObjValEDNtupleProducer::~ObjValEDNtupleProducer>:" << std::endl;
  std::cout << " numEvents processed = " << numEvents_processed_ << std::endl;

  std::cout << "done." << std::endl;
}

void ObjValEDNtupleProducer::beginJob() {}

void ObjValEDNtupleProducer::produce(edm::Event& evt, const edm::EventSetup& es)
{
  //std::cout << "<ObjValEDNtupleProducer::produce>:" << std::endl;

//--- add run, luminosity section and event number
  std::auto_ptr<edm::RunNumber_t> runNumberPtr(new edm::RunNumber_t(evt.id().run()));
  evt.put(runNumberPtr, "run");
  std::auto_ptr<edm::LuminosityBlockNumber_t> lumiSectionPtr(new edm::LuminosityBlockNumber_t(evt.id().luminosityBlock()));
  evt.put(lumiSectionPtr, "lumisection");
  std::auto_ptr<edm::EventNumber_t> eventNumberPtr(new edm::EventNumber_t(evt.id().event()));
  evt.put(eventNumberPtr, "event");

  typedef std::vector<double> vdouble;

//--- compute values to be filled in ntuple
  for ( std::vector<ntupleEntryType*>::iterator ntupleEntry = ntupleEntries_.begin();
	ntupleEntry != ntupleEntries_.end(); ++ntupleEntry ) {
    if ( !(*ntupleEntry)->objValExtractor_ ) {
      edm::LogError ("ObjValEDNtupleProducer::fillntuplees")
	<< " No ObjValExtractor set for ntuple = " << (*ntupleEntry)->ntupleName_
	<< " --> skipping !!";
      continue;
    }
    std::auto_ptr<double> toPut = std::auto_ptr<double>(new double);
    *toPut = (*(*ntupleEntry)->objValExtractor_)(evt);
    evt.put(toPut, (*ntupleEntry)->ntupleName_);
  }
  
  for ( std::vector<ntupleVectorEntryType*>::iterator ntupleEntry = ntupleVectorEntries_.begin();
	ntupleEntry != ntupleVectorEntries_.end(); ++ntupleEntry) {
    if ( !(*ntupleEntry)->objValExtractor_ ) {
      edm::LogError ("ObjValEDNtupleProducer::fillntuplees")
	<< " No ObjValExtractor set for ntuple = " << (*ntupleEntry)->ntupleName_
	<< " --> skipping !!";
      continue;
    }
    
    vdouble values = (*(*ntupleEntry)->objValExtractor_)(evt);

    if ( (*ntupleEntry)->indices_.size() == 0 ) {
      try { 
	std::auto_ptr<vdouble> toPut = std::auto_ptr<vdouble>(new vdouble((*(*ntupleEntry)->objValExtractor_)(evt)));
	evt.put(toPut, (*ntupleEntry)->ntupleName_);
      } catch ( cms::Exception e ) { 
	edm::LogError("ObjValEDNtupleProducer::produce")
	  << " ObjValExtractor plugin name = " << (*ntupleEntry)->ntupleName_  << " caused exception --> rethrowing !!";
	throw e;
      }
    } else {
      try { 
	vdouble values = (*(*ntupleEntry)->objValExtractor_)(evt);
      } catch ( cms::Exception e ) { 
	edm::LogError("ObjValEDNtupleProducer::produce")
	  << " ObjVectorValExtractor plugin name = " << (*ntupleEntry)->ntupleName_  << " caused exception --> rethrowing !!";
	throw e;
      }	
      //std::cout << " values = " << format_vdouble(values) << std::endl;
      std::auto_ptr<vdouble> toPut = std::auto_ptr<vdouble>(new vdouble());
      for ( vunsigned::const_iterator index = (*ntupleEntry)->indices_.begin(); 
	    index != (*ntupleEntry)->indices_.end(); ++index ) {
	//std::cout << "index = " << (*index) << std::endl;
	if ( (*index) < values.size() ) toPut->push_back(values[*index]);
      }
      //std::cout << "--> toPut = " << format_vdouble(*toPut) << std::endl;
      evt.put(toPut, (*ntupleEntry)->ntupleName_);
    }
  }

  ++numEvents_processed_;
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(ObjValEDNtupleProducer);
