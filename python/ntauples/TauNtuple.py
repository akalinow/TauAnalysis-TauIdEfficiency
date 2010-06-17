import string
import hashlib

class TauNtuple(object):
    " Class that maps user-friendly variable names to complicated TTree branches "
    def __init__(self, ntuple_name, collection_name, variables=[]):
        self.collection = collection_name
        self.variables = list(variables)
        # Build
        self.mapping_dict = dict(
            (var, "%(ntuple_name)s#%(collection_name)s#%(var)s" % {
                 'ntuple_name':ntuple_name, 
                 'collection_name':collection_name, 
                 'var':var } 
            ) for var in variables)

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
        return TauNtupleExpression(self.substitute(expression)) 

class TauNtupleExpression(object):
    def __init__(self, value):
        self.value = value
    def hash_string(self):
        " Return a hash string corresponding to this expression " 
        return hashlib.md5(self.value).hexdigest()
    def __hash__(self):
        " Return a hash integer corresponding to this expression " 
        return int(self.hash_string(), 16)
    def check_other(self, other):
        # To be implemented if needed later
        pass
        #if isinstance(other, TauNtupleExpression):
        #    if self.events is not other.events:
        #       raise ValueError
    def __add__(self, other):
        self.check_other(other)
        return TauNtupleExpression(
            "(%s) + (%s)" % (self.value, other))
    def __sub__(self, other):
        self.check_other(other)
        return TauNtupleExpression(
            "(%s) - (%s)" % (self.value, other))
    def __mul__(self, other):
        self.check_other(other)
        return TauNtupleExpression("(%s)*(%s)" % (self.value, other))
    def __div__(self, other):
        self.check_other(other)
        return TauNtupleExpression("(%s)/(%s)" % (self.value, other))
    def __pow__(self, other):
        if not isinstance(other, int) or other < 2:
            raise ValueError
        output = self
        for i in range(other-1):
            output = output*output
        return output
    def __and__(self, other):
        self.check_other(other)
        return TauNtupleExpression("(%s) && (%s)" % (self.value, other))
    def __or__(self, other):
        self.check_other(other)
        return TauNtupleExpression("(%s) || (%s)" % (self.value, other))
    def __lt__(self, other):
        self.check_other(other)
        return TauNtupleExpression("(%s) < (%s)" % (self.value, other))
    def __le__(self, other):
        self.check_other(other)
        return TauNtupleExpression("(%s) <= (%s)" % (self.value, other))
    def __gt__(self, other):
        self.check_other(other)
        return TauNtupleExpression("(%s) > (%s)" % (self.value, other))
    def __ge__(self, other):
        self.check_other(other)
        return TauNtupleExpression("(%s) >= (%s)" % (self.value, other))
    def __eq__(self, other):
        self.check_other(other)
        return TauNtupleExpression("(%s) == (%s)" % (self.value, other))
    def __str__(self):
        return self.value
    def false(self):
        " Negate expression "
        return TauNtupleExpression(
             "!(%s)" % self.value)

if __name__ == "__main__":
    # Some tests
    # Some mutable object 
    pt = TauNtupleExpression('pt')
    eta = TauNtupleExpression('eta')
    new = pt | eta
    print new

