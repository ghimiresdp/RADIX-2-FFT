# list of all imports
import math
import cmath
import matplotlib.pyplot as plt


# Class Radix2
class Radix2:
    # initialization of all the required variables
    def __init__(self, n):

        self.N = n  # Size of an array
        self.stages = int(math.log(self.N, 2))  # no. of stages/iterations in butterfly structure
        self.xn = [0 + 0j] * self.N     # initialization of an array with default values 0+0j
        self.xn1 = self.xn.copy()  # to store bit reversed values
        self.XN = self.xn.copy()    # to store the final FFT value

        # print(self.N, ' Points DFT Initialized !!')
        # print('no. of stages = ', self.stages)

        # Twiddle Factor
        self.w = cmath.exp((-1j) * 2 * math.pi / self.N)
        self.W = []
        for i in range(int(self.N / 2)):
            self.W.append(self.w ** i)

    # Reversing bits for DIT in First stage.
    def reverse_bits(self, n):
        result = 0
        for i in range(self.stages):
            result <<= 1
            result |= n & 1
            n >>= 1
        return result

    # assign the input values
    def inputs(self, input_list):
        for i in range(input_list.__len__()):
            self.xn[i] = complex(input_list[i])
            self.xn1[self.reverse_bits(i)] = self.xn[i]

    # Print out FFT values
    def print_values(self):
        print('\ninput x[n]:')
        for i in range(self.xn.__len__()):
            print('x[', i, ']:\t', self.xn[i])
        print('\noutput X(N):')
        for i in range(self.XN.__len__()):
            print('X(', i, '):\t', self.XN[i])

        print('\n Twiddle Factors W(N)^n:')
        for i in range(int(self.N / 2)):
            print('value of w(', self.N, ')^', i, ': ', self.W[i])

    # Perform the FFT Radix-2 Algorithm
    def run(self):
        for i in range(1, self.stages + 1):     # iteration through stages
            xn_tmp = [0 + 0j] * self.N  # assign a temporary value for calculation

            # Iteration throughout the length of input array
            # in which selection depends upon the stage
            for j in range(0, self.N, 2 ** i):
                index = []  # Selected Indices of the butterfly for operation

                for k in range(j, j + 2 ** i):
                    index.append(k)
                # print(index)
                mid_point = index.__len__() // 2
                power = 0
                step = math.log(self.N, 2) - (i - 1)

                # perform FFT operation in the selected indices
                # for example:
                #   -> a+W^n*b
                #   -> a-W^n*b
                # to understand what the values might be,
                # please have a look at the butterfly structure
                for m in range(mid_point):
                    xn_tmp[index[m]] = \
                        self.xn1[index[m]] + (self.w ** power) * self.xn1[index[m + mid_point]]

                    xn_tmp[index[m + mid_point]] = \
                        self.xn1[index[m]] - (self.w ** power) * self.xn1[index[m + mid_point]]

                    # print(index[m], xn_tmp[index[m]])
                    # print(index[m + mid_point], xn_tmp[index[m + mid_point]])
                    power += step

            # Replace the current output to previous input for next stage operation
            self.xn1 = xn_tmp.copy()

        # Copy the final output for in the variable XN
        self.XN = self.xn1.copy()

    def result(self):
        return self.XN

    # Plot the real and imaginary values in a MatPlotLib graph
    # please have a look at the Python MatPlotLib first
    # if you have not installed MatPlotLib, please install using:
    # pip install matplotlib
    def plot(self, input_values=False, output_values=False, grid=True):
        xn_real = []
        xn_img = []
        XN_real = []
        XN_img = []

        # Separate the real and imaginary values to plot in the graph
        for i in range(self.N):
            xn_real.append(self.xn[i].real)
            xn_img.append(self.xn[i].imag)
            XN_real.append(self.XN[i].real)
            XN_img.append(self.XN[i].imag)

        # For both the Input and the Output Waveforms
        if input_values is True and output_values is True:
            plt.subplot(2, 1, 1)
            plt.plot(range(self.N), xn_real, color='blue')
            plt.plot(range(self.N), xn_img, color='red')
            plt.xlabel('[n]')
            plt.ylabel('x[n]')
            if grid is True:
                plt.grid()

            plt.subplot(2, 1, 2)
            plt.plot(range(self.N), XN_real, color='blue')
            plt.plot(range(self.N), XN_img, color='red')
            plt.xlabel('[n]')
            plt.ylabel('X(N)')
            if grid is True:
                plt.grid()

        # For both the Input Waveform only
        elif input_values is True:
            plt.plot(range(self.N), xn_real, color='blue')
            plt.plot(range(self.N), xn_img, color='red')
            plt.xlabel('[n]')
            plt.ylabel('x[n]')
            if grid is True:
                plt.grid()

        # For the Output Waveform only
        elif output_values is True:
            plt.plot(range(self.N), XN_real, color='blue')
            plt.plot(range(self.N), XN_img, color='red')
            plt.xlabel('[n]')
            plt.ylabel('X(N)')
            if grid is True:
                plt.grid()
        plt.show()

