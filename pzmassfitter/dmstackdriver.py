#!/usr/bin/env python
####################
# Example driver for maxlike masses in Clusters project
####################

import maxlike_controller, maxlike_bentstep_voigt as mbv
import  astropytable_filehandler as atf, maxlike_masses


#######################

makeController = lambda : maxlike_controller.Controller(modelbuilder =  mbv.BentVoigtShapedistro(),
                                    filehandler = atf.AstropyTableFilehandler(),
                                    runmethod = maxlike_masses.SampleModelToFile())


controller = makeController()



if __name__ == '__main__':

    controller.run_all()

