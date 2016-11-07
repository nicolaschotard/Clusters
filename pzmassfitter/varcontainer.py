###############
# Basic Python Object that is a dual between an object and a dictionary
################

__cvs_id__ = "$Id$"

################


class VarContainer(dict):

    def __getattr__(self, name):

        return self[name]

    def __setattr__(self, name, value):

        self[name] = value

