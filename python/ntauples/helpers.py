import hashlib
import ROOT

'''

Helper functions for the ntaupls package.

Author: Evan K. Friis (UC Davis)

'''

def make_mean_rms_pave(plot, x_low=0.6, y_low=0.85, x_high=0.99, y_high=0.99):
    mean = plot.GetMean()
    rms = plot.GetRMS()
    output = ROOT.TPaveText(x_low, y_low, x_high, y_high, "brNDC")

def copy_aliases_from(from_tree, to_tree):
   aliases = [ alias.GetName() for alias in from_tree.GetListOfAliases() ]
   for alias in aliases:
      to_tree.SetAlias(alias, from_tree.GetAlias(alias))

def copy_aliases(tchain):
   ''' Ensure that all TTrees in a TChain have their aliases set correctly '''
   # Get the first ttree in the chain
   tchain.LoadTree(0)
   temp_tree = tchain.GetTree()
   # Copy the aliases
   copy_aliases_from(temp_tree, tchain)

def filter_aliases(aliases, *match_to):
    " Yields aliases from aliases whose beginning fileds match "
    for alias in aliases:
        #fields = alias.split('#')
        matches = True
        for field, match_it in zip(alias, match_to):
            if field != match_it:
                matches = False
                break
        if matches:
            yield alias

def make_unique_name(*args):
    " Make a unique hash using the concatenation of str(args) "
    return hashlib.md5(''.join(str(arg) for arg in args)).hexdigest()
    #return "adfasdfsdafasdf"



