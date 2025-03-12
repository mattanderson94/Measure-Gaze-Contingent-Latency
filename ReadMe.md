# A method for measuring closed-loop latency in gaze-contingent rendering without extra equipment.

## 📌 Description
This is a method for precisely estimating the closed-loop latency of gaze-contingent rendering. By closed-loop latency, we mean the time between a change in where an observer is looking on a stimulus, and a corresponding update to the stimulus. This is critical for artifact-free foveated rendering, and smooth experiments. Various methods exist, but ours has the advantage of requiring no additional hardware beyond what is typically used in eye tracking studies (a monitor and an eye tracker).

Full details are in our paper [here (to come)](https://google.com)

## 🚀 Installation

```bash
# Clone the repository
git clone https://github.com/mattanderson94/Measure-Gaze-Contingent-Latency.git

```
Or simply download a zip and extract to a local folder.

This MATLAB code is written for an Eyelink 1000+ and has only been tested on a display computer configured as follows:
```MATLAB
% -----------------------------------------------------------------------------------------------------
% MATLAB Version: 23.2.0.2668659 (R2023b) Update 9
% MATLAB License Number: xxxxxxx
% Operating System: Microsoft Windows 10 Pro Version 10.0 (Build 19045)
% Java Version: Java 1.8.0_202-b08 with Oracle Corporation Java HotSpot(TM) 64-Bit Server VM mixed mode
% -----------------------------------------------------------------------------------------------------
% MATLAB                                                Version 23.2        (R2023b)    
% Psychtoolbox                                          Version 3.0.19      17 February
```

This code should be compatible with some other versions of MATLAB and Psychtoolbox, but this is not guaranteed. For proper Psychtoolbox installation instructions if needed, [click here](http://psychtoolbox.org/download).

You may also need to download the latest [Eyelink development kit](https://www.sr-research.com/support/forum-9.html).

This code also requires a correctly installed version of the [Edf2Mat© Matlab Toolbox](https://github.com/uzh/edf-converter). Follow installation instructions on their github page. In a few of the scripts, you will be prompted to add the ./@Edf2Mat/ folder to the current MATLAB search path, so make sure the path is added before running any analysis scripts. 

## 📖 Usage

To measure gaze contingent latency, run 
```MATLAB
RunPTBClosedLoopTimingTest();
```

To analyse the data to extract the latency estimates, run
```MATLAB
EstimateLatencyFromEyeTrackerData();
```

And if you want to validate the method, you can compare expected versus measured latencies by running
```MATLAB
CompareExpectedVersusObservedLatencies();
```
The logic of this validation approach is described in detail in [our paper](https://google.com)

## 🤝 Contributing

Contributions are welcome! Please follow these steps:
1. Fork the repository
2. Create a new branch (`git checkout -b feature-branch`)
3. Commit your changes (`git commit -m 'Add some feature'`)
4. Push to the branch (`git push origin feature-branch`)
5. Open a Pull Request


## 📞 Contact
For any questions, feel free to reach out:
- Email: matt.anderson@berkeley.edu

