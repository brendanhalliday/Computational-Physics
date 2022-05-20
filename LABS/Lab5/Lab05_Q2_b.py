# This program was written by Brendan Halliday
# For PHY407 lab 5 question 2 part b.

# This program uses fourier transform and applies a low pass
# filter on a wav audio file.
# This is done by setting all sound frequencies above 880Hz
# to zero and leaving in all the lower frequencies.
# It then returns the filtered audio file as a bass
# bosted version of the original. 


from scipy.io.wavfile import read, write
import numpy as np
import matplotlib.pyplot as plt
import os # not one of the sanctioned libraries, however it is necessary

# Sometimes reading the wav file does not work.
# So, I've added these two lines to find directory of
# the running python file and then change the
# current working directory to the directory of the 
# running file.
# This is because sometimes the current working
# directory is different from the one in which the file 
# is running.

directory_name = os.path.dirname(__file__) # find directory of current running python file
os.chdir(directory_name) # change the working directory to the directory that includes the running file

# Now, provided that the .wav file is included in the same directory as this python file,
# all should run smoothly.


#----------------------------------------------------------------|



# read the data into two stereo channels
# sample is the sampling rate, data is the data in each channel,
# dimensions [2, nsamples]
sample, data = read(r"GraviteaTime.wav")
# sample is the sampling frequency, 44100 Hz
# separate into channels
channel_0 = data[:, 0]
channel_1 = data[:, 1]


dt = 1/sample  # s

N_points = len(channel_0)
# length of interval
T = N_points*dt
# convert to (angular) frequency
w_freq = np.arange(N_points/2+1)* 2 * np.pi / T
print(len(w_freq))
# each sample occurs ever 1/sample seconds
time = dt * np.arange(0, N_points)


# plot channel_0
plt.subplot(2,1,1)
plt.plot(time, channel_0)
plt.ylabel('Channel_0 signal')
plt.title("Channel_0 and Channel_1 signals vs. time", y =1.08)
plt.grid()
plt.xlim(1.000, 1.005)
plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))

# plot channel_1
plt.subplot(2,1,2)
plt.plot(time, channel_1)
plt.xlabel('time (s)')
plt.ylabel('Channel_1 signal')
plt.grid()
plt.xlim(1.000, 1.005)
plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
plt.savefig("rawsound.png", dpi = 300)
plt.show()

# Perform a FOurier transform on both signals
A0 = np.fft.rfft(channel_0)
AMPLITUDE0 = np.abs(A0)
A1 = np.fft.rfft(channel_0)
AMPLITUDE1 = np.abs(A1)


# plot raw frequency data
plt.subplot(2,1,1)
plt.plot(w_freq, AMPLITUDE0, label="Channel 0")
plt.ylabel('|Fourier Coefficients|')
plt.title('Raw Frequency Spectrum - Fourier Transform of audio signal', y=1.08)
plt.xlim(0, 40000)
plt.grid()
plt.legend()

plt.subplot(2,1,2)
plt.plot(w_freq, AMPLITUDE1,label="Channel 1")
plt.xlabel('Angular Frequency (Hz)')
plt.ylabel('|Fourier Coefficients|')
plt.xlim(0, 40000)
plt.grid()
plt.legend()
plt.savefig("rawfreq.png", dpi =300)
plt.show()

A0[abs(w_freq) > 880.0] = 0.
A1[abs(w_freq) > 880.0] = 0.

# This plot plots the new filtered frequency data
plt.subplot(2,1,1)
plt.plot(w_freq, abs(A0), label="Channel 0")
plt.ylabel('|Fourier Coefficients|')
plt.title("Filtered Frequency Spectrum - Fourier Transform of audio signal", y = 1.08)
plt.xlim(0, 1600)
plt.legend()
plt.grid()

plt.subplot(2,1,2)
plt.plot(w_freq, abs(A1), label="Channel 1")
plt.xlabel('Angular Frequency (Hz)')
plt.ylabel('|Fourier Coefficients|')
plt.xlim(0, 1600)
plt.grid()
plt.legend()
plt.savefig("filteredfreq.png", dpi = 300)
plt.show()


Back0 = np.fft.irfft(A0)
Back1 = np.fft.irfft(A1)

plt.subplot(2,1,1)
plt.plot(time, Back0)
plt.ylabel('Channel 0 Audio Signal')
plt.title('Filtered Audio Signal', y = 1.08)
plt.grid()
plt.xlim(1.,1.005)
plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))

plt.subplot(2,1,2)
plt.plot(time, Back1)
plt.xlabel('time (s)')
plt.ylabel('Channel 1 Audio Signal')
plt.grid()
plt.xlim(1.,1.005)
plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
plt.savefig("filteresaudio.png", dpi = 300)
plt.show()


data_out = np.empty(data.shape, dtype = data.dtype)
# fill data_out
data_out[:, 0] = Back0 
data_out[:, 1] = Back1
write('GraviteaTime_lpf.wav', sample, data_out)