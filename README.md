# TS_Dopamine_Optogenetics
Dopamine sensor and D1 or D2 neuron calcium signals in the tail of striatum with optogenetic manipulation of dopamine.

We recorded from dopamine axons in the tail of striatum (TS) using a dopamine sensor DA3m and calcium activity from D1 or D2 neurons medium spiny neurons in TS using a calcium indicator GCaMP8m . Recordings were performed when head-fixed mice experience a sensory stimulus, optogenetic activation of TS dopamine axons, or both at the same time.

**File structure 
The datasets are divided into 5 experiments with different task conditions. Each file contains the data for all mice used in one experiment.
‘0_masking' or '0_no_masking’  : Activity of TS dopamine axons was recorded with a dopamine sensor DA3m when mice experienced optogenetic activation of dopamine axons with or without a background masking light. The masking light was used to eliminate the neuron activity caused by optogenetics excitation light itself. Two red lamps were placed on both sides of the mouse with approximately 15 cm distance to their eyes, and were turned on before the mice entered the arena and kept on throughout the experiments. All the following experiments were performed with the masking light on. 
‘1_natural_novel_response’: Activity of TS dopamine axons was recorded with a dopamine sensor DA3m when mice were presented with a novel multi-modal sensory stimulus (blue LED light and a 100db complex tone ). Each sensory pulse lasted for 0.5 s and was repeated four times with 0.5 s inter-stimulus intervals.
‘2_opto_three_frequencies’: Activity of TS dopamine axons was recorded with a dopamine sensor DA3m when dopamine axons were optogenetically activated at different frequencies (10, 20 and 40Hz).
‘3_task_xx_response_dxbx’: Dopamine sensor (DA3m) signals in TS or calcium sensor (GCaMP8m) signals in TS D1 or D2 neurons were recorded while mice were repeatedly presented with the same sensory stimulus ( blue LED light and a 50db complex tone ) with or without dopamine axon optogenetic activation. D1 and D2 activities were recorded for two days with two blocks per day. Block 1 includes 20 trials with no optogenetic stimulation (“no stim trials”), and block 2 includes 30 trials, 15 trials with optogenetic stimulation (“stim trials”) and 15 “no stim-trials” in a pseudo-random order. The blocks are indicated by “dxbx” in the file name. For example, “d1b2” refers to day 1, block 2  .
‘4_opto_only_D1_response’: Activity of TS D1 neurons was recorded with a calcium sensor GCaMP8m when dopamine axons were optogenetically activated at 20Hz.

**Data structure in each file
In each .mat file, there is a 1xn cell array called "signal_raw", where n is the number of mice.
For experiments 0, 1, 2, the first half are control mice without opsin and the second half are mice with opsin. The order of mice is the same across files.
For experiment 3, the first half are mice with opsin and the second half are control mice without opsin. The order of mice is the same within D1 recordings and within D2 recordings.
For experiment 4, all mice are with opsin.

In the cell array, each element is a tx3 matrix that contains the recorded data of one mouse from one session. t is the length of that session, and the columns are three separate channels recorded in unit of Volt. The sampling rate is 1000Hz.
The first column is the recorded calcium activity  , which is converted to electric signal by a photodetector.
The second column is the signal that controls optogenetic stimulation laser.
The third column is the signal that controls the sensory stimulus (blue LED light and tone).

