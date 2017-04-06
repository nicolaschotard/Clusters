#!/usr/bin/env python
####################
# Example driver for maxlike masses in Clusters project
####################

import maxlike_controller
import maxlike_bentstep_voigt
import maxlike_masses
import astropytable_filehandler as atf


#######################

makeTestingController = lambda: maxlike_controller.Controller(modelbuilder=maxlike_masses.LensingModel(),
                                                              filehandler=atf.AstropyTableFilehandler(),
                                                              runmethod=maxlike_masses.ScanModelToFile())

# Use definition below [BentVoigtShapedistro()] if no STEP2 shear calibration required 
makeController = lambda: maxlike_controller.Controller(modelbuilder=maxlike_bentstep_voigt.BentVoigtShapedistro(),
                                                       filehandler=atf.AstropyTableFilehandler(),
                                                       runmethod=maxlike_masses.SampleModelToFile())

# Use definition below [BentVoigt3Shapedistro()] to enable STEP2 shear calibration
# makeController = lambda: maxlike_controller.Controller(modelbuilder=maxlike_bentstep_voigt.BentVoigt3Shapedistro(),
#                                                      filehandler=atf.AstropyTableFilehandler(),
#                                                       runmethod=maxlike_masses.SampleModelToFile())


controller = makeController()


if __name__ == '__main__':

    controller.run_all()

