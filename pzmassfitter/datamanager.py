#!/usr/bin/env python
######################
# Class to manage data needed to build the model
# and to track cuts made to said files
# only general enough to make sense to maxlike_masses.py
######################

from __future__ import with_statement
import os
import pwd
import time
import datetime
import inspect
import getpass
from . import util

#######################


class ManagedDataException(Exception):
    pass

#######################


def pretty_localtime(t): return time.asctime(time.localtime(t))


def pretty_datetime(t): return time.asctime(t.timetuple())

#########


class HistoryCollection(list):

    ###

    def comment(self, comment):

        self.append(CommentHistoryItem(comment))

    ###

    def log(self):
        seperator = '*****************\n'
        return seperator + seperator.join([x.log() for x in self]) + seperator

    ###

    def dump(self, file):

        with open(file, 'w') as output:
            output.write(self.log())

    ###

    def __str__(self):
        return self.log()

#########


class HistoryItem(object):

    def __init__(self, action):

        self.action = action
        self.time = datetime.datetime.now()

    ###

    def __str__(self):

        return self.log()

    ###

    def log(self):

        header = '%(action)s  at  %(time)s\n' \
            % {'action': self.action, 'time': pretty_datetime(self.time)}

        body = self._log()

        return header + body

#####


class InitHistoryItem(HistoryItem):

    def __init__(self):

        super(InitHistoryItem, self).__init__('INIT')

    ###

    def _log(self):

        str = '''
DataManager Initialized by %(user)s
''' % {'user': getpass.getuser()}

        return str

######


class StoreHistoryItem(HistoryItem):

    def __init__(self, name, replace=False):
        action = 'STORE'
        if replace is True:
            action = 'REPLACE'
        super(StoreHistoryItem, self).__init__(action)

        self.name = name

    ###

    def _log(self):

        if isinstance(self.name, str):
            return '\n%s initialized.\n\n' % self.name

        return '\n' + '\n'.join(['%s initialized.' % x for x in self.name]) + '\n\n'

#######


class StoreFunctionHistoryItem(HistoryItem):

    def __init__(self, name, method, methodargs, methodkeywords, replace=False):

        action = 'STORE'
        if replace is True:
            action = 'REPLACE'

        super(StoreFunctionHistoryItem, self).__init__(action)

        self.name = name
        self.method = method
        self.methodargs = methodargs
        self.methodkeywords = methodkeywords
        self.methodsrcfile = inspect.getfile(method)
        self.methodsrcstat = os.stat(self.methodsrcfile)
        self.methodsrcode = inspect.getsource(method)

    ####

    def _log(self):

        if isinstance(self.name, str):
            printedName = self.name
        else:
            printedName = ', '.join(self.name)

        return '''
%(name)s initialzed
    Using function %(method)s 
     Called with args = %(args)s
        and keywords = %(keywords)s
     Defined in %(srcfile)s   Modifed on %(srcmoddate)s
     Looks like:
%(srccode)s
''' % {'name': printedName,
            'method': self.method.__name__,
            'args': str(self.methodargs),
            'keywords': str(self.methodkeywords),
            'srcfile': self.methodsrcfile,
            'srcmoddate': pretty_localtime(self.methodsrcstat.st_mtime),
            'srccode': self.methodsrcode}

######


class OpenHistoryItem(HistoryItem):

    def __init__(self, name, file, method, methodargs, methodkeywords):
        super(OpenHistoryItem, self).__init__('OPEN')

        self.name = name
        self.file = file
        self.method = method
        self.methodargs = methodargs
        self.methodkeywords = methodkeywords
        self.filestat = os.stat(file)
        self.methodsrcfile = inspect.getfile(method)
        self.methodsrcstat = os.stat(self.methodsrcfile)

    #######

    def _log(self):

        str = '''
%(name)s <==  %(file)s
     Owned by %(owner)s   Modified on %(moddate)s
    Using function %(method)s 
     Called with args = %(args)s
      and keywords = %(keywords)s
     Defined in %(srcfile)s   Modifed on %(srcmoddate)s
''' % {'name': self.name,
            'file': os.path.abspath(self.file),
            'owner': pwd.getpwuid(self.filestat.st_uid)[0],
            'moddate': pretty_localtime(self.filestat.st_mtime),
            'method': self.method.__name__,
            'args': self.methodargs,
            'keywords': self.methodkeywords,
            'srcfile': self.methodsrcfile,
            'srcmoddate': pretty_localtime(self.methodsrcstat.st_mtime)}

        return str

#######


class UpdateHistoryItem(HistoryItem):

    def __init__(self, name, method, args, keywords):
        super(UpdateHistoryItem, self).__init__('UPDATE')

        self.name = name
        self.method = method
        self.methodargs = args
        self.methodkeywords = keywords
        self.methodsrcfile = inspect.getfile(method)
        self.methodsrcstat = os.stat(self.methodsrcfile)
        self.methodsrcode = inspect.getsource(method)

    ####

    def _log(self):

        if isinstance(self.name, str):
            printedName = self.name
        else:
            printedName = ', '.join(self.name)

        return '''
%(name)s modified
    Using function %(method)s 
     Called with args = %(args)s
      and keywords = %(keywords)s
     Defined in %(srcfile)s   Modifed on %(srcmoddate)s
     Looks like:
%(srccode)s
''' % {'name': printedName,
            'method': self.method.__name__,
            'args': self.methodargs,
            'keywords': self.methodkeywords,
            'srcfile': self.methodsrcfile,
            'srcmoddate': pretty_localtime(self.methodsrcstat.st_mtime),
            'srccode': self.methodsrcode}

#####################


class CommentHistoryItem(HistoryItem):

    def __init__(self, comment):
        super(CommentHistoryItem, self).__init__('COMMENT')

        self.comment = comment

    ####

    def _log(self):

        return '\n%s\n\n' % self.comment


#######################

class DataManager(util.VarContainer):

    def __init__(self):

        super(DataManager, self).__init__(self)

        self.__dict__['history'] = HistoryCollection()
        self.history.append(InitHistoryItem())

    ############################

    def __setitem__(self, name, value):

        self._addItem(name, value)

    #############################

    def __setattr__(self, name, value):

        self.history.append(StoreHistoryItem(name))
        self._addItem(name, value)

    #############################

    def _addItem(self, name, value, replace=False):

        if name not in self or replace is True:
            super(DataManager, self).__setitem__(name, value)
        else:
            raise ManagedDataException

    #############################

    def comment(self, comment):

        self.history.append(CommentHistoryItem(comment))

    #############################

    def open(self, name, file, open, *args, **keywords):

        self.history.append(OpenHistoryItem(name=name,
                                            file=file,
                                            method=open,
                                            methodargs=args,
                                            methodkeywords=keywords))

        self._addItem(name, open(file, *args, **keywords))

    ##############################

    def store(self, name, src, *args, **keywords):

        self._store(name, src, replace=False, args=args, keywords=keywords)

    ##############################

    def _store(self, name, src, replace=False, args=None, keywords=None):

        if args is None:
            args = ()
        if keywords is None:
            keywords = {}

        try:
            tostore = src(*args, **keywords)
            self.history.append(StoreFunctionHistoryItem(name=name,
                                                         method=src,
                                                         methodargs=args,
                                                         methodkeywords=keywords,
                                                         replace=replace))

        except Exception as e:
            tostore = src
            self.history.append(StoreHistoryItem(name, replace=replace))

        if isinstance(name, str):
            self._addItem(name, tostore, replace=replace)

        else:
            for curname, curobj in zip(name, tostore):
                self._addItem(curname, curobj, replace=replace)

    ###############################

    def replace(self, name, src, *args, **keywords):

        self._store(name, src, replace=True, args=args, keywords=keywords)

    ###############################

    def update(self, filterFunction, toUpdate, *args, **keywords):

        if len(args) == 0 and len(keywords) == 0:
            try:
                theFilter = filterFunction(self)
            except TypeError:
                theFilter = filterFunction()
        else:
            theFilter = filterFunction(*args, **keywords)

        self.history.append(UpdateHistoryItem(name=toUpdate.keys(),
                                              method=filterFunction,
                                              args=args,
                                              keywords=keywords))

        for item, updateFunc in toUpdate.iteritems():

            self._addItem(item, getattr(self[item], updateFunc)(
                theFilter), replace=True)

    ######################################

    def listappend(self, name, othermanager):

        if name not in self:
            super(DataManager, self).__setitem__(name, [])

        self.__getattr__(name).append(othermanager)

        self.history.extend(othermanager.history[1:])
