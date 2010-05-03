def copy_aliases(tchain):
   ''' Ensure that all TTrees in a TChain have their aliases set correctly '''
   # Get the first ttree in the chain
   tchain.LoadTree(0)
   temp_tree = tchain.GetTree()
   # Copy the aliases
   aliases = [ alias.GetName() for alias in temp_tree.GetListOfAliases() ]
   for alias in aliases:
      tchain.SetAlias(alias, temp_tree.GetAlias(alias))

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

