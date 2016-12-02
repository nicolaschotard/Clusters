#!/usr/bin/env python
######################
# @file bashreader.py
# @author Douglas Applegate
# @date 2/26/08
#
# @brief Interprets simple bash scripts to parse them for variables
#     This way python can share existing config files.
#
# This is meant as a library. Test routines will run if the script is executed.
########################

from __future__ import with_statement
import unittest
import re
import datetime
import os


class DoesNotParseException(Exception):
    pass

######################################################

class BashConfig(dict):
    '''Reads a simple bash file and returns a dictionary-interface
         object containing the bash variables set in the script.
       It will also process any bash export statements.
    '''

    def __init__(self, bashstr=None):
        self.parse(bashstr)

    ###################

    def __getattr__(self, attr):
        try:
            return self[attr]
        except KeyError:
            raise AttributeError('Unknown Variable: %s' % attr)

    ###################

    def _parseInt(self, bashstr):

        if re.match(r'^(\d+)$', bashstr):
            return int(bashstr)

        raise DoesNotParseException

    ###################

    def _parseFloat(self, bashstr):

        if re.match(r'^[-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?$', bashstr):
            return float(bashstr)

        raise DoesNotParseException

    ###################

    _variableRE = r'\${?(\w+)(?:\[(\d+)\])?}?'

    def _parseVar(self, bashstr):

        match = re.match(r'^%s$' % BashConfig._variableRE, bashstr)
        if match is None:
            raise DoesNotParseException

        return self._replaceVar(match)

    ###################

    def _replaceVar(self, match):

        varname = match.group(1)
        if str.lower(varname) in self:
            val = getattr(self, str.lower(varname))
        elif varname in os.environ:
            val = os.environ[varname]
        else:
            raise NameError

        index = match.group(2)
        if index is None:
            return val

        return val[int(index)]

    ###

    def _parseString(self, bashstr):

        ###

        def subVars(match):
            return str(self._replaceVar(match))

        ###

        strippedQuotes = bashstr.strip('"\'')
        return re.sub(BashConfig._variableRE, subVars,
                      strippedQuotes)

    ##################

    def _parseTime(self, bashstr):

        match = re.match(r'(\d+):(\d+):(\d+)', bashstr)
        if match is not None:
            hours = int(match.group(1))
            min = int(match.group(2))
            sec = int(match.group(3))
        
            return datetime.time(hours, min, sec)

        raise DoesNotParseException

    ###################

    def _parseArray(self, bashstr):

        match = re.match(r'\((.+)\)', bashstr)
        if match is None:
            raise DoesNotParseException

        arrayDef = match.group(1)
        return [self._parseAtomic(x) for x in arrayDef.split()]

    ###################

    def _parseDict(self, bashstr):

        entries = re.findall(r'\[(\d+)\]\s*=\s*(.+?)[\s\)]', bashstr)

        if len(entries) == 0:
            raise DoesNotParseException

        return dict([(int(x), self._parseAtomic(y) ) for x, y in entries])

    ##################

    def _parseEval(self, bashstr):

        match = re.match(r'\$\(\((.+)\)\)', bashstr)
        if match is None:
            raise DoesNotParseException

        toEval = self._parseString(match.group(1))

        return eval(toEval)

    ##################

    def _parseAtomic(self, bashstr):

        return self._runParsers(bashstr, parsers=[self._parseInt,
                                                  self._parseFloat,
                                                  self._parseVar,
                                                  self._parseString])

    ##################

    def _runParsers(self, bashstr, parsers):

        for parser in parsers:
            try:
                return parser(bashstr)
            except DoesNotParseException:
                continue

    ##################

    def _parseExport(self, bashstr):

        match = re.match(r'export ((\w+)(=.+)?)', bashstr)
        if match is None:
            raise DoesNotParseException

        if match.group(3) is not None:
            self._parseAssignment(match.group(1))

        attr = match.group(2)
        os.environ[attr] = str(getattr(self, str.lower(attr)))

    ##################

    def _parseAssignment(self, bashstr):

        match = re.match(r'(\w+)(?:\[(\w+)\])?=(.+)', bashstr)
        if match is None:
            raise DoesNotParseException

        attr = str.lower(str.strip(match.group(1)))
        rawIndex = match.group(2)
        rawVal = str.strip(match.group(3))

        valParsers = [self._parseInt,
                      self._parseFloat,
                      self._parseVar,
                      self._parseTime,
                      self._parseEval,
                      self._parseDict,
                      self._parseArray,
                      self._parseString]

        val = self._runParsers(rawVal, valParsers)

        if rawIndex is None:
            self[attr] = val
        else:
            index = int(rawIndex)
            if attr in self:
                self[attr][index] = val
            else:
                self[attr] = {index: val}

    ###################

    def _parseReadFile(self, bashstr):

        match = re.match(r'\. (.+)', bashstr)
        if match is None:
            raise DoesNotParseException

        filename = match.group(1)
        self.parseFile(filename)

    ###################

    def _parseSemicolon(self, bashstr):

        substatements = bashstr.split(';')
        if len(substatements) == 1:
            raise DoesNotParseException

        for substatement in substatements:
            self._parseLine(substatement)

    ###################

    def _parseLine(self, bashstr):

        lineParsers = [self._parseSemicolon,
                       self._parseReadFile,
                       self._parseAssignment,
                       self._parseExport]

        self._runParsers(bashstr.strip(), lineParsers)

    ###################

    def parse(self, bashstr):

        if bashstr is None:
            return

        lines = map(str.strip, bashstr.splitlines())

        for line in lines:

            self._parseLine(line)

    ###################

    def parseFile(self, filename):

        with open(filename) as input:
            for line in input:
                self._parseLine(line)

#######################################################
#######################################################

######################
# USER METHODS
##############

def parse(bashstr):

    result = BashConfig(bashstr)
    return result

#######################################################

def parseFile(filename):

    result = BashConfig()
    result.parseFile(filename)
    return result

######################################################
######################################################

############################
# TESTING CLASSES
##################

class TestParseBash(unittest.TestCase):

    def testReadIntVar(self):

        config = parse('NFRAMES=20')
        self.assertEquals(config.nframes, 20)
        self.assertEquals(type(config.nframes), type(5))

    #############

    def testReadFloat(self):

        config = parse('OBSLAT=19.82861111')
        self.assertAlmostEquals(config.obslat, 19.82861111, 8)

    #############

    def testReadStr(self):

        config = parse('INSTRUMENT=SUBARU')
        self.assertEquals(config.instrument, 'SUBARU')

    #############

    def testReadTime(self):

        config = parse('REFERENCETIME=22:00:00')
        self.assertEquals(config.referencetime, datetime.time(22, 0, 0))

    ##############

    def testReadIntIntMap(self):

        config = parse('OVSCANX1=([6]=1  [7]=1  [3]=1  [4]=1  [9]=1  [8]=2055)')
        self.assertEquals(config.ovscanx1, {6:1, 7:1, 3:1, 4:1, 9:1, 8:2055})

    ##############

    def testReadMultiplelines(self):

        bash = '''OBSLAT=19.82861111
OBSLONG=155.48055556

REFERENCETIME=22:00:00
'''
        config = parse(bash)
        self.assertAlmostEquals(config.obslat, 19.82861111, 8)
        self.assertAlmostEquals(config.obslong, 155.48055556, 8)
        self.assertEquals(config.referencetime, datetime.time(22, 0, 0))

    ##############

    def testReadReference(self):

        bash = '''REF=5

TEST=$REF
TEST2=${REF}
'''
        config = parse(bash)
        self.assertEquals(config.ref, 5)
        self.assertEquals(config.test, 5)
        self.assertEquals(config.test2, 5)

    ##############

    def testReadReferenceInConcat(self):

        bash = '''BIN=/home/bin
P_READLINK=${BIN}/readlink'''

        config = parse(bash)
        self.assertEquals(config.bin, '/home/bin')
        self.assertEquals(config.p_readlink, '/home/bin/readlink')

    ###############

    def testBadReference(self):

        self.assertRaises(NameError, parse, 'TEST=$BAD')

    ##############

    def testEnvReference(self):

        config = parse('TEST=$HOME')
        self.assertEquals(config.test, os.environ['HOME'])

    ##############

    def testMapRef(self):

        bash = '''STATSALLIM=([1]=1000 [2]=2000 [3]=1000 [4]=1000)

TEST=$STATSALLIM[1]
TEST2=${STATSALLIM[1]}'''

        config = parse(bash)
        self.assertEquals(config.statsallim, {1:1000, 2:2000, 3:1000, 4:1000})
        self.assertEquals(config.test, 1000)
        self.assertEquals(config.test2, 1000)

    ##############

    def testNonIntMap(self):

        bash = 'TEST=([1]="a" [2]="b")'

        config = parse(bash)
        self.assertEquals(config.test, {1:'a', 2:'b'})

    ##############

    def testStripQuotes(self):

        config = parse('TEST="blab"')
        self.assertEquals(config.test, 'blab')

    ##############

    def testParseArray(self):

        config = parse('TEST=(1 2 3 4)')
        self.assertEquals(config.test, [1, 2, 3, 4])

    ##############

    def testParseEvaluate(self):

        bash = '''STATSALLIM=([1]=1000 [2]=2000 [3]=1000 [4]=1000)

STATSXMIN=$(( ${STATSALLIM[1]} - ${STATSALLIM[3]} / 2 ))
STATSXMAX=$(( ${STATSALLIM[1]} + ${STATSALLIM[3]} / 2 ))
STATSYMIN=$(( ${STATSALLIM[2]} - ${STATSALLIM[4]} / 2 ))
STATSYMAX=$(( ${STATSALLIM[2]} + ${STATSALLIM[4]} / 2 ))

'''
        config = parse(bash)
        self.assertEquals(config.statsallim, {1:1000, 2:2000, 3:1000, 4:1000})
        self.assertEquals(config.statsxmin, 500)
        self.assertEquals(config.statsxmax, 1500)
        self.assertEquals(config.statsymin, 1500)
        self.assertEquals(config.statsymax, 2500)

    ##################

    def testExportVar(self):

        bash = '''TEST=5
export TEST'''
        config = parse(bash)
        self.assertTrue('TEST' in os.environ)
        self.assertEquals(os.environ['TEST'], '5')

        del os.environ['TEST']

    #################

    def testExportVar2(self):
        config = parse('export TEST=13')
        self.assertEquals(config.test, 13)
        self.assertTrue('TEST' in os.environ)
        self.assertEquals(os.environ['TEST'], '13')

        del os.environ['TEST']

    ##################

    def testExportIssue1(self):

        config = parse('export TEMPDIR=.')
        self.assertEquals(config.tempdir, '.')
        self.assertTrue('TEMPDIR' in os.environ)
        self.assertEquals(os.environ['TEMPDIR'], '.')



    #################

    def testReadFile(self):

        cfile = 'test.ini'

        if not os.path.exists(cfile):
            output = open(cfile, 'w')
            output.write('TEST=5\n')
            output.close()

        config = parse('. test.ini')
        self.assertEquals(config.test, 5)

        os.remove(cfile)

    ##################

    def testSemicolons(self):

        config = parse('TEST=5 ; export TEST; TEST2=10')
        self.assertEquals(config.test, 5)
        self.assertTrue('TEST' in os.environ)
        self.assertEquals(os.environ['TEST'], '5')
        self.assertEquals(config.test2, 10)

    ###################

    def testSetDefinedMapElement(self):

        bash = '''TEST=([1]=5 [2]=10)
                  TEST[1]=8'''
        config = parse(bash)
        self.assertEquals(config.test, {1:8, 2:10})

    ####################

    def testSetUndefinedMapElement(self):

        config = parse('TEST[1]=3')
        self.assertEquals(config.test, {1:3})

    ####################

    def testParseFile(self):

        cfile = 'test.ini'

        if not os.path.exists(cfile):
            output = open(cfile, 'w')
            output.write('TEST=5\n')
            output.close()

        config = parseFile(cfile)
        self.assertEquals(config.test, 5)

        os.remove(cfile)

    ###########################


######################################################

def test():

    testcases = [TestParseBash]
    suite = unittest.TestSuite(map(unittest.TestLoader().loadTestsFromTestCase,
                                   testcases))
    unittest.TextTestRunner(verbosity=2).run(suite)

##############################################################

if __name__ == '__main__':

    test()

