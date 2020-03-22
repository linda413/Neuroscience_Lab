%% Preprocessing
%     * Before you start, ensure all scripts and data files are on your Matlab path!
%     * Put NCS files on ThetaCruncher 
%         * IT IS IMPORTANT THAT FILES ARE CONSISTENTLY NAMED AND PROPERLY ORGANIZED BEFORE CONVERSION (the next step).
%         * This way the files all get named appropriately and put in appropriately named folders.
%         * A little data organization and planning goes a really long way in the end. 
%         * All of the following analyses use the animal name and task/recording session that is listed in the file name.
%           So, a good rule of naming for this lab seems to be: animalName_recordingSession_dateCollected
%         * Note: For those recording sessions that have multiple animals contained in one session (and whose NSC files 
%                 are labeled CSC 1:64), it might be best to separate those into individual animal folders before conversion. 
%                 You can use changeNames.m to do this quickly, but it requires heavy modification. So, be careful!
%     * Use convertNCS2MATBatch_AS.m to convert NCS files to MAT files.
%         * Follow instructions in script and enter UDIs
%         * Calls convertNCS2MAT.m and Nlx2MatCSC (code provided by Neuralynx)
%     * Create a Bad Channels sheet in your Animal Codes excel file (see template or prior study's file)
%         * The subsequent steps use this "bad channels" sheeet to load only those files that have good data, 
%           which cuts down on computation time and ensures you don't include bad data in your analyses.
%         * To detect bad channels you can either:
%             * 1) Run seeAllChannels.m to quickly detect which channels are noise
%                 * Follow instructions in script and enter UDIs
%             * 2) Skip the above step and proceed to the next step of artifact rejection, where you can also detect bad channels.
%         * Keep track of those bad channels in your excel file
%     * Run readAnimalCodes.m to ensure you are reading in your Animal Information excel file correctly
%     * Artifact Rejection
%         * For this step, artifact rejection is mainly for large artifact that indicates a seizure or a tether disconnection.
%           It should be present across all recording channels and usually reaching the amplitdue cutoff of 3 mV.
%           For epileptic animals, do not mark the spikes; they will be detected and quantified in the next step.
%         * Run one of the scripts below.
%         * Run DAVIS_MultiTrace_Manual_Artifact_Reject.m (script by Matt at UofA)
%             * Enter in command line: DAVIS_MultiTrace_Manual_Artifact_Reject([],450,1)
%                 * If you leave the first argument empty (i.e. []), then use the GUI to select the CSC files for that animal and task.
%                 * The second argument assigns the size of the time window to view the data.
%         * Run artifactRejectionBatch_AS.m
%             * Enter in the UDIs and the script then calls DAVIS_MultiTrace_Manual_Artifact_Reject.m
%             * It just cuts down on time by bypassing the GUI and looping through the animals.
%         * Both scripts load the good channels for each animal for a recording session, 
%           and the user marks the global artifact present across all channels. 
%         * The marked intervals are stored in the original data file and structure but in a field called bad_intervals. 
%           Later on, other scripts will use this information to remove power estimates for those bad segments.
%
%% Spike Counting
%     * This section is only important if you are working with epileptic animals not TBI.
%     * Run countInterictalSpikesBatch_AS.m
%         * Calls countInterictalSpikes_AS.m
%           * This function has a few UDIs at the beginning that you could/should look at/modify!
%             * Suggested parameters: numSTDsignal = 4.5, absoluteMin = 0.5, numQuantile = 0.9
%           * Calls waveletDecomp_AS.m
%     * Run plotNORspikeTotals.m to view the results from the automated spike counting
%     * Run periSpikeHistogram.m to get an average spike shape for each animal
%
%% Behavior and Performance
%   * Run plotNORbehavior.m
%       * This code is highly specific to the NOR task and the first two minutes of behavior.
%   * Run plotNORperformance.m
%       * Calls readPerformance.m
%       * It also relies on the spike counts from the previous step.
%   * Run plotSpikesAllTimes.m
%       * This code plots the spike counts over multiple recording sessions. It is quite experiment/study specific.
%
%% Theta Analyses
%     * Power
%         * Use powerAnalyses.m
%             * Calls readAnimalCodes.m, plotNORspikeTotals.m, readBadChannels.m, power_AS.m (waveletDecomp_AS.m, addMirroredBuffers.m, interpolateBadData.m), readNORbehavior.m, plotBoxScatter.m
%             * Outputs multiple plots
%     * Pepisode
%         * To calculate and save the output from pEpisode, run pEpisodeBatch_AS.m. Then, pEpisodeAnalyses.m can uses these saved output files to average across groups and create plots. These saved files just help on computation time.
%         * Use pepisodeAnalyses.m
%             * Calls readAnimalCodes.m, plotNORspikeTotals.m, readBadChannels.m, pepisode_AS.m (waveletDecomp_AS.m, addMirroredBuffers.m, interpolateBadData.m), readNORbehavior.m, plotBoxScatter.m
%             * Outputs multiple plots
%     * Coherence
%         * Use coherenceAnalyses.m
%             * Calls  coherence_AS.m, mscohereNaN.m
%     * Compare NOR1 and NOR2
%         * Use interactingAnalyses.m
%             * Calls readAnimalCodes.m, plotNORspikeTotals.m, readBadChannels.m
%     * Phase
%         * Use phaseAnalyses.m (still being developed)
%
%% Extra Code
%     * spikeCluster.m: Allows you to sort the animals by spike rate rather than group/treatment
%     * changeName.m: Allows you to rename and move multiple animal recording session files (like for the long term recordings) into single files. This code needs to be heavily modified, so please use with caution.
%     * removeExtraEEG.m: Some recordings don?t immediately start at ?time point 0,? so you may want to remove/truncate the beginning or end of the recording to match up with your actual behavior.
%
%% Notes
%     * There is a button/command called "Publish" that allows you to view code and code outputs as a PDF or HTML file.