import string

class TauNtuple(object):
    " Class that maps user-friendly variable names to complicated TTree branches "
    def __init__(self, events, ntuple_name, collection_name, variables=[]):
        self.events = events
        self.collection = collection_name
        self.variables = list(variables)
        # Build
        self.mapping_dict = dict(
            (var, 
             "%(ntuple_name)s#%(collection_name)s#%(var)s" % {
                 'ntuple_name':ntuple_name, 
                 'collection_name':collection_name, 
                 'var':var} )
            for var in variables)

    def substitute(self, expr):
        " Retrive the correct TTree draw()-type string for a given variable "
        return string.Template(expr).substitute(self.mapping_dict)

    def __repr__(self):
        " Print helpful information"
        output = "Tau Ntuple - collection:%s\n" % self.collection
        for var in self.variables:
            output += " * * %s\n" % var
        return output

    def expr(self, expression):
        " Build an expression from this ntuple "
        return TauNtupleExpression(self.events, self.substitute(expression)) 

class TauNtupleExpression(object):
    def __init__(self, events, value):
        self.events = events
        self.value = value
    def check_other(self, other):
        if isinstance(other, TauNtupleExpression):
            if self.events is not other.events:
                raise ValueError
    def __add__(self, other):
        self.check_other(other)
        return TauNtupleExpression(
            self.events, 
            "(%s) + (%s)" % (self.value, other))
    def __sub__(self, other):
        self.check_other(other)
        return TauNtupleExpression(
            self.events, "(%s) - (%s)" % (self.value, other))
    def __mul__(self, other):
        self.check_other(other)
        return TauNtupleExpression(
            self.events, "(%s)*(%s)" % (self.value, other))
    def __div__(self, other):
        self.check_other(other)
        return TauNtupleExpression(
            self.events, "(%s)/(%s)" % (self.value, other))
    def __and__(self, other):
        self.check_other(other)
        return TauNtupleExpression(
            self.events, "(%s) && (%s)" % (self.value, other))
    def __or__(self, other):
        self.check_other(other)
        return TauNtupleExpression(
            self.events, "(%s) || (%s)" % (self.value, other))
    def __lt__(self, other):
        self.check_other(other)
        return TauNtupleExpression(
            self.events, "(%s) < (%s)" % (self.value, other))
    def __le__(self, other):
        self.check_other(other)
        return TauNtupleExpression(
            self.events, "(%s) <= (%s)" % (self.value, other))
    def __gt__(self, other):
        self.check_other(other)
        return TauNtupleExpression(
            self.events, "(%s) > (%s)" % (self.value, other))
    def __ge__(self, other):
        self.check_other(other)
        return TauNtupleExpression(
            self.events, "(%s) >= (%s)" % (self.value, other))
    def __eq__(self, other):
        self.check_other(other)
        return TauNtupleExpression(
            self.events, "(%s) == (%s)" % (self.value, other))
    def __str__(self):
        return self.value
    def false(self):
        " Negate expression "
        return TauNtupleExpression(
            self.events, "!(%s)" % self.value)

if __name__ == "__main__":
    # Some tests
    # Some mutable object 
    events = ['a', 'b', 'c']
    pt = TauNtupleExpression(events, 'pt')
    eta = TauNtupleExpression(events, 'eta')
    new = pt | eta
    print new

