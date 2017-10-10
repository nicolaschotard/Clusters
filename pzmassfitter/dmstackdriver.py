"""Example driver for maxlike masses in Clusters project."""


from . import maxlike_controller as mcont
from . import maxlike_bentstep_voigt
from . import maxlike_masses
from . import astropytable_filehandler as atf


def makeTestingController():
    return mcont.Controller(modelbuilder=maxlike_masses.LensingModel(),
                            filehandler=atf.AstropyTableFilehandler(),
                            runmethod=maxlike_masses.ScanModelToFile())


# Use definition below [BentVoigtShapedistro()] if no STEP2 shear calibration required
def makeController():
    return mcont.Controller(modelbuilder=maxlike_bentstep_voigt.BentVoigtShapedistro(),
                            filehandler=atf.AstropyTableFilehandler(),
                            runmethod=maxlike_masses.SampleModelToFile())


# Use definition below [BentVoigt3Shapedistro()] to enable STEP2 shear calibration
# Obsolete: now flag in BentVoigtShapedistro() to use or not STEP2 shear calibration
#def makeController_sc():
#    return mcont.Controller(modelbuilder=maxlike_bentstep_voigt.BentVoigt3Shapedistro(),
#                            filehandler=atf.AstropyTableFilehandler(),
#                            runmethod=maxlike_masses.SampleModelToFile())


controller = makeController()
#controller = makeController_sc()


if __name__ == '__main__':

    controller.run_all()
