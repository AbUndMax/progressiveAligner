# Progressive Aligner for Multiple Sequence Alignment (MSA)

## Features
* Consensus sequence pairwise alignment (adapted Needleman-Wunsch) as profile-profile alignment technique 
* uses either Neighbour Joining to calculate a guiding tree which determines the order in which two profiles get aligned.
* or newly calculated distances (using consensus sequences) between each profile to determine which profiles to align next.

This program is used as Batch with the following command:

first compile the progressiveAlignment.java file
```Bash
javac ProgressiveAlignment.java
```

start the program with the following parameters:
```Bash
java ProgressiveAlignment <Match-Score> <Mis-Match-Score> <Gap-Penalty> <Path-To-FastaFile>
```

an optional verbose mode is available:
```Bash
<-v>
```

For the current version of the program only consensus Profile-Profile comparison is implemented!

### The output could look like this:
```
#### The MSA is computed by comparing !CONSENSUS! sequences!

MKRKLISEQSLKELNSTPGHEITEQFHRTK
M--KLISEQSLKELNSTPGHEKFGEKRVLD
----LISEQ--IEIKIST--PGHEFHR-TK
----LISEQ--IEIKI--------FHR-TK
----LISEQ-HFDLN---FH------R-TK
 .. *****   .             . ..
```