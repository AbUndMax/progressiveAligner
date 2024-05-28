# Progressive Aligner for Multiple Sequence Alignment (MSA)

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