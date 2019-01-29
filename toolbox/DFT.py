import numpy
import pyopencl as cl
import reikna.cluda as cluda
from reikna.fft import FFT as fftcl



def DFT(inputsignal):
    signal = numpy.array(inputsignal, dtype = numpy.complex128)
    device = cl.get_platforms()[1].get_devices(device_type=cl.device_type.GPU)
    ctx = cl.Context(device)
    queue = cl.CommandQueue(ctx)
    
    api = cluda.ocl_api()
    thr = api.Thread(queue)
    signal_dev = thr.to_device(signal)
    fft_res_dev = thr.array(signal.shape, dtype = numpy.complex128)
    
    FFT = fftcl(signal_dev).compile(thr)
    FFT(fft_res_dev, signal_dev)
    
    res = fft_res_dev.get()
    return res
	
def IDFT(inputsignal):
    
    signal = numpy.array(inputsignal, dtype = numpy.complex128)
    device = cl.get_platforms()[1].get_devices(device_type=cl.device_type.GPU)
    ctx = cl.Context(device)
    queue = cl.CommandQueue(ctx)
    
    api = cluda.ocl_api()
    thr = api.Thread(queue)
    
    signal_dev = thr.to_device(signal)
    ifft_res_dev = thr.array(signal.shape, dtype = numpy.complex128)
    
    IFFT = fftcl(ifft_res_dev).compile(thr)
    IFFT(ifft_res_dev, signal_dev, inverse = True)
    
    res = ifft_res_dev.get()
    return res