# jeevan-thanu-icarus-challenge
APPROACH:-

STEP1:I started to implement the bpsk modulation and demodulation using input bitstream as a random sequence of 0s and 1s.I then generated the code for it using chatgpt.

STEP2:I researched on chatgpt the proper order of implementation of all phases mentioned in the pdf.Then I proceeded according to the order.

STEP3:According to the order first I gave the phase2 dataset for snr and told chatgpt to write the code for snr calibration and then combine the bpsk code and the snr code.Then I gave the phase 1 dataset for timing asked it to integrate the timing code with the rest.Then I gave the phase 4(doppler) and phase 3 (coding) dataset respectively and asked it to integrate the entire thing to one code.

STEP4:Then I matched with the challenge table and it showed that the Reed-solomon values are missing.So I added the value as RS(15,11)..and asked it to modify the entire code.

STEP5:Then I told it to write the code for idle python.Then I downloaded numpy,scipy,panda and imported thyem.Then I copy pasted the code in idle python and then ran it.Then after debugging many errors related to location of files and folder creation;the excel sheet was getting generated.

STEP6:Then when I was matching the excel data produced to the performance threshold given in the challenge table ;it was not at all matching.Then after many modifications to the code I jumped back to the code produced after after step 4 and started working on it again.

STEP7:Then debugged each phase again and combined it.Then when I matched with challenge table ..it was perfectly matching.And all the rz.npy,meta.json of the samples of all phases were getting decoded.

STEP8:Then I generated the ber vs snr graph through chatgpt.Then I generated the code for constellation diagrams and doppler diagrams and then I combined the original phase1-4 code the constellation diagram code and the doppler diagram code.

STEP9:Then I modified the code so that it saves the constellation diagram plots,doppler diagram plots in a folder named plot directly in the desktop.

Key challenges and lessons learnt:-

1.First I was getting an error of no value of RS PASSED .So I gave the value as RS(15,11) and it was working.And even in the challenge table it was specified that we had to take RS values as 15,11.

2.When I copy and pasted the code in python(phase1-4 code)..and ran it .A empty excel file was getting generated. In the ICARUS folder.Then after debugging and many prompts I realised that the all the phase folders must be extracted in the location where  the excel file is getting created which is the ICARUS folder.Then after I ran the code still it was not working.Then I realised that the elements like 0db,5db,10db,15db were directly getting extracted into the ICARUS folder.Then I created folder for each phase ..and gave it the same name as in the code.Then I put all the respective elements of the phases inside their respective folders.And then the excel sheet with values was getting generated.

3.When I was generating codes for phase 1-4 and executing on idle python.It was not at all running after a couple of times .And it showed an permission denied error.Then I realised that after you execute the code one time and when you open the excel sheet ,u must close it cause windows cant replace or change the new values in the excel sheet if it is open.

4.The values generated in the excel sheet was not matching with the challenge table provided.After many prompts I realised that I have to go back to the original code after I added RS(15,11) as during that run phase 1 and 2 were working.Then after woring on that phase 3and 4 were working but 1 and 2 had errors.I was always getting naïve result and not the performance threshold.Then I debugged each and every phase individually and after I tuned the code by many prompts and then when I compared it to the challenge table ,it was working and the values were matching with the performance threshold.

5.Then I used chatgpt to plot the ber vs snr graph.Then I realised that I could plot the ber vs snr graph through a code in the main code itself.Then after many tries I realised that I couldn’t get accurate values for the grapgh through the python code as it was unable to take all the realtime data which was getting created in the excel sheet and plot the graph.The graph shape was also not the same.Therefore I decided to go ahead by feeding the entire excel sheet to chatgpt and come up with the ber vs snr graph.As it was more accurate.


TOOLS AND RESOURCES
CHATGPT
https://chatgpt.com/share/68b9d8ab-e704-800b-b799-db27e9f94fcb
