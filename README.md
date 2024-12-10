[![GitHub](https://img.shields.io/badge/GitHub-progressiveAligner-b07219?logo=github)](https://github.com/AbUndMax/progressiveAligner)
[![Java](https://img.shields.io/badge/Java-11+-b07219)](https://openjdk.org/projects/jdk/11/)
[![License](https://img.shields.io/badge/License-CC_BY--NC_4.0-blue)](https://github.com/AbUndMax/progressiveAligner/blob/main/LICENSE.md)
[![Badge](https://img.shields.io/github/v/release/AbUndMax/progressiveAligner?color=brightgreen)](https://github.com/AbUndMax/progressiveAligner/releases/latest)

# Progressive Aligner for Multiple Sequence Alignment (MSA)

## Features
* Consensus sequence pairwise alignment (adapted Needleman-Wunsch) as profile-profile alignment technique
* Uses either Neighbour Joining to calculate a guiding tree which determines the order in which two profiles get aligned.
* Or newly calculated distances (using consensus sequences) between each profile to determine which profiles to align next.

## Run the application:
Download the JAR file from the [latest release](https://github.com/AbUndMax/progressiveAligner/releases/latest).
Navigate to the folder in which the JAR was saved in and run it with:
```Bash
java -jar progressiveAligner.jar <arguments...>
```

## Available Parameters

| **Parameter**          | **Short Flag** | **Type**      | **Required** | **Description**                                                                                                                                             | **Default** |
|-------------------------|----------------|---------------|--------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------|-------------|
| `--fastaPath`           | `-fp`         | `[s] String`  | Mandatory    | Specifies the path to the FASTA file containing at least 2 sequences.                                                                                       | -           |
| `--matchScore`          | `-ms`         | `[i] Integer` | Optional     | Positive value of the match score.                                                                                                                          | 4           |
| `--misMatchScore`       | `-mms`        | `[i] Integer` | Optional     | Positive value of the mismatch score.                                                                                                                       | 2           |
| `--gapPenalty`          | `-gp`         | `[i] Integer` | Optional     | Positive value of the gap penalty.                                                                                                                          | 1           |

## Available Commands

| **Command**          | **Short Flag** | **Description**                                                                                                                                       |
|-----------------------|----------------|-------------------------------------------------------------------------------------------------------------------------------------------------------|
| `Consensus`          | `c`            | Uses newly computed distances between profiles, based on consensus sequences, to decide which profiles to align next.                                  |
| `NeighbourJoining`   | `nj`           | Builds a guiding tree using the Neighbour Joining method to determine the order of profile-profile alignments.                                         |

## Clone and work on the SourceCode:

### 1. Prerequisites
Ensure the following tools are installed:
- **Java Development Kit (JDK)**: Version 11 or higher (21 recommended).
- **Maven**: To build and manage dependencies.
- **Git**: To clone the repository.

### 2. Clone the Repository

### 3. Build the Project
Use Maven to compile the code, run tests, and package the application:
```Bash
mvn clean install
```
This will:
- Resolve dependencies.
- Compile the source code.
- Run unit tests (if applicable).
- Create a packaged JAR file in the target directory.