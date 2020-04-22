from info_bslop0_apr4highsnr import basic, paths, config, zbins, plotting, source_nofz_pars, sysmaps
from xipm_xcorrmod_binslop0_apr4highsnr import GGL, Measurement, ResponsesScale, ResponsesProjection

T = True
F = False

run_measurement = T


if run_measurement:
    print 'Starting measurement class...'
    #gglensing = GGL(basic, config, paths)
    measurement = Measurement(basic, config, paths, zbins, plotting)
    if not basic['plot_blinded']:
        measurement.run()
        #measurement.save_boostfactors_2pointfile() #deprecated, without errors
        #measurement.save_2pointfile('gt')


