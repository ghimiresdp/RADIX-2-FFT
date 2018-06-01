# This is a sample program usage for the fft library
# To understand the details about the radix-2 FFT, please look at the class fft

# import the class Radix2
import fft

# define a  as 8 point Input
a = fft.Radix2(64)

# define input values
# value = [1 / pow(2, 0.5),
#          1,
#          1 / pow(2, 0.5),
#          0,
#          -1 / pow(2, 0.5), -1, -1 / pow(2, 0.5), -0]
value = [1]*63
a.inputs(value)

# perform FFT Operations
a.run()
# print(a.result())
a.print_values()
a.plot(input_values=True, output_values=True, grid=True)
