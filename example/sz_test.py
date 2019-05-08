import time
import pysz
import numpy as np

data = np.random.rand(1000,1000)
bound = 1e-4
compressor = pysz.Compressor((pysz.ConfigBuilder()
    .errorBoundMode(pysz.ABS)
    .absErrBound(1e-4)
    .build()))

t1 = time.process_time()
compressed = compressor.Compress(data)
t2 = time.process_time()
decompressed = compressor.Decompress(compressed, data.shape, data.dtype)
t3 = time.process_time()

print("compression time", t2-t1)
print("decompression time", t3-t2)
