#!/usr/bin/env python
####################
# Example driver for maxlike masses in Clusters project
####################

import maxlike_controller, maxlike_bentstep_voigt
import  astropytable_filehandler as atf, maxlike_masses


#######################

makeTestingController = lambda : maxlike_controller.Controller(modelbuilder =  maxlike_masses.LensingModel(),
                                    filehandler = atf.AstropyTableFilehandler(),
                                    runmethod = maxlike_masses.ScanModelToFile())

makeController = lambda : maxlike_controller.Controller(modelbuilder =  maxlike_bentstep_voigt.BentVoigtShapedistro(),
                                    filehandler = atf.AstropyTableFilehandler(),
                                    runmethod = maxlike_masses.SampleModelToFile())


controller = makeController()



if __name__ == '__main__':

    controller.run_all()

