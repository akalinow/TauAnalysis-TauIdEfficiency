#include "TauAnalysis/TauIdEfficiency/plugins/ObjValEDNtupleProducer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

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

    //typedef std::vector<std::string> vstring;
    //vstring cfgInputTagNames = cfgObjValExtractor.getParameterNamesForType<edm::InputTag>();

    // Determine if we want the whole collection or just specific indices.  Default case, index=0
    typedef std::vector<unsigned> vunsigned;
    vunsigned indices; 
    
    if ( cfgObjValExtractor.exists("indices") ) {
      vunsigned cfgIndices = cfgObjValExtractor.getParameter<vunsigned>("indices");
      for ( vunsigned::const_iterator index = cfgIndices.begin(); index != cfgIndices.end(); ++index ) {
	indices.push_back(*index);
      }
    } else if ( cfgObjValExtractor.exists("vector") && cfgObjValExtractor.getParameter<bool>("vector") == true ) {
      // if we want all objects, indicate it by leaving indices empty
    } else {
      // default case = take first index
      indices.push_back(0);
    }
    
    // Get the names of all the columns
    edm::ParameterSet columnsCfg = cfgObjValExtractor.getParameter<edm::ParameterSet>("columns");
    std::vector<std::string> columnNames = columnsCfg.getParameterNamesForType<std::string>();
    
    // Build each column source
    for( std::vector<std::string>::const_iterator columnName = columnNames.begin();
	 columnName != columnNames.end(); ++columnName ) {
      std::string columnCfg = columnsCfg.getParameter<std::string>(*columnName);
      
      // concatentate relavant information into a single PSet to send to the ObjValExtractors
      edm::ParameterSet tempNtupleCfg(cfgObjValExtractor);
      //edm::ParameterSet tempNtupleCfg;
      //tempNtupleCfg.addParameter<edm::InputTag>("src", cfgCollection);
      //for ( vstring::const_iterator cfgInputTagName = cfgInputTagNames.begin();
      //      cfgInputTagName != cfgInputTagNames.end(); ++cfgInputTagName ) {
      //  edm::InputTag cfgInputTag = cfgObjValExtractor.getParameter<edm::InputTag>(*cfgInputTagName);
      //  tempNtupleCfg.addParameter<edm::InputTag>(*cfgInputTagName, cfgInputTag);
      //}
      tempNtupleCfg.addParameter<std::string>("value", columnCfg);
      
      // Check if we are taking specific indices
      if ( indices.size() ) {
	// Loop over each index
	for ( vunsigned::const_iterator index = indices.begin(); index != indices.end(); ++index ) {
	  edm::ParameterSet tempNtupleCfgWithIndex = tempNtupleCfg;
	  tempNtupleCfgWithIndex.addParameter<unsigned>("index", *index);
	  // Build entry
	  ObjValExtractorBase* objValExtractor = 
	    ObjValExtractorPluginFactory::get()->create(pluginTypeObjValExtractor, tempNtupleCfgWithIndex);
	  std::string niceColumnName;
	  if ( indices.size() == 1 && (*index) == 0 )
	    niceColumnName = createBranchName(ntupleName_, cfgCollection.label(), *columnName);
	  else
	    niceColumnName = createBranchName(ntupleName_, cfgCollection.label(), *columnName, *index);
	  ntupleEntryType* ntupleEntry = new ntupleEntryType(niceColumnName, objValExtractor); 
	  ntupleEntries_.push_back(ntupleEntry);
	}
      } else {
	// Build a vector getter
	ObjValVectorExtractorBase* objValExtractor = 
	  ObjValVectorExtractorPluginFactory::get()->create(pluginTypeObjValExtractor, tempNtupleCfg);
	std::string niceColumnName = createBranchName(ntupleName_, cfgCollection.label(), *columnName);
	ntupleVectorEntryType* ntupleEntry = new ntupleVectorEntryType(niceColumnName, objValExtractor);
	ntupleVectorEntries_.push_back(ntupleEntry);
      }
    }
  }
  
  // Register all the products with the framework
  for ( std::vector<ntupleEntryType*>::const_iterator entry = ntupleEntries_.begin();
	entry != ntupleEntries_.end(); ++entry ) {
    std::string name = (*entry)->ntupleName_;
    produces<double>(name).setBranchAlias(name);
  }
  
  for ( std::vector<ntupleVectorEntryType*>::const_iterator entry = ntupleVectorEntries_.begin();
	entry != ntupleVectorEntries_.end(); ++entry ) {
    std::string name = (*entry)->ntupleName_;
    produces<std::vector<double> >(name).setBranchAlias(name);
  }
  
  numEvents_filled_ = 0;
  
  std::cout << "done." << std::endl;
}

ObjValEDNtupleProducer::~ObjValEDNtupleProducer()
{
  std::cout << "<ObjValEDNtupleProducer::~ObjValEDNtupleProducer>:" << std::endl;
  std::cout << " numEvents filled into EDMNtuples = " << numEvents_filled_ << std::endl;
  
  std::cout << "done." << std::endl;
}

void ObjValEDNtupleProducer::beginJob() {}

void ObjValEDNtupleProducer::produce(edm::Event& evt, const edm::EventSetup& es)
{
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

    std::auto_ptr<vdouble> toPut = std::auto_ptr<vdouble>(new vdouble((*(*ntupleEntry)->objValExtractor_)(evt)));

    evt.put(toPut, (*ntupleEntry)->ntupleName_);
  }
  
  ++numEvents_filled_;
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(ObjValEDNtupleProducer);
