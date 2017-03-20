# ProductionNetworkSimulator
This code can simulate how shocks on firms are propagted through a supply chain network. The code is based on the model proposed by [Hallegatte, 2008; Henriet and Hallegatte, 2008]. The model is modified to simulate actual production networks by [Inoue and Todo, 2017].

* [Hallegatte, 2008] Hallegatte, S. (2008). An adaptive regional input-output model and its application to the assessment of the economic cost of katrina. Risk analysis, 28(3):779-799.
* [Henriet and Hallegatte, 2008] Henriet, F. and Hallegatte, S. (2008). Assessing the consequences of natural disasters on production networks: a disaggregated approach.
* [Inoue and Todo, 2017] Inoue, H. and Todo, Y. (2017) Propagation of Negative Shocks Through Firm Networks: Evidence from Simulation on Comprehensive Supply-Chain Data. Available at SSRN: https://ssrn.com/abstract=2932559

# 1. How to start using this?

Most readers who use this model is probably not familiar with compile or gitHub. If you are windows 64bit user, please follow the steps below. If not, you need to complie the source files.

For windows 64bit user:
1. Find and click "release" link on this page. Or click this. https://github.com/HiroyasuInoue/ProductionNetworkSimulator/releases
2. Find "Downloads" list and download proNet.exe and Source code (either of zip or tar is OK).
3. Uncompress source code. They are source codes and sample data. Put them in the same directory with proNet.exe
4. Open windows command line and execute proNet.exe.
5. A simplest execution example is

   proNet.exe toyATableLine.txt toyCVector.txt toyPiniVector.txt toyFirmAffiliation.txt outputPactToy.txt

   or try

   proNet.exe -d toyDelta01.txt -f outputFCToy.txt -s outputSTATtoy.txt toyATableLine.txt toyCVector.txt toyPiniVector.txt toyFirmAffiliation.txt outputPactToy.txt

6. Check Usage section for more information.

# 2. Usage

## 2.1 Minimum usage

Following is the minimum usage of proNet.exe

proNet.exe aTableLineFile cVectorFile piniVectorFile firmAffiFile pactOutputFile

The five args are input and output files and they are the minimum fileset. All files should not have any contradiction. If there is a contradiction between files, the results are not guaranteed completely. Formats for the five files. See toyX files of source code.

1. aTableLineFile (input)

(supplier firm ID) (client firm ID) (amount (normally money))

This shows the amount of supply chain network.

2. cVectorFile (input)

(supplier firm ID) (amount (normally money))

This shows the amount of supply toward final consumer by a firm.

3. piniVectorFile (input)

(supplier firm ID) (amount (normally money))

This shows the amount of supply initially provided by a firm.

4. firmAffiFile (input)

(firm ID) (affiliation ID)

This shows affiliations of firms.

5. pactOutputFile

(time) (total production of firms)

This shows total production of firms at each time step.

## 2.2 Usage with options

proNet.exe has a bunch of options. They are put between proNet.exe and aTableLineFile. The format is

proNet.exe [option] aTableLineFile cVectorFile piniVectorFile firmAffiFile pactOutputFile

[option]s are following.
* -a (int): adaptation (0: as much as possible, 1: only once)
* -A (str): adaptation counter file (available only when -a is on)
* -d (str): delta file (default null)
* -f (str): final consumption output file (default null)
* -F : fulldebug output
* -o (int): order algorithm type. 0: normal, 1: keep initial demand, 2: ignore negative inv adjustment
* -p (file:int:int...): get snapshot of production at indicated step(s). for more than one snapshot, int should be separated by :
* -r (int): random seed (default 331)
* -R (int): rationing algorithm 0. proportional, 1. FC low priority, 2. lower has high priority, 3. 1 and 2 (default proportional)
* -s (str): stats output file (default null)
* -t (int): simulation step (default 10)
* -v (str): value added output file (default null)

The detail follows as below.

* -a (int): adaptation (0: as much as possible, 1: only once)

(under implementation)

* -A (str): adaptation counter file (available only when -a is on)

(under implementation)

* -d (str): delta file (default null)

Delta file can specifies delta for firm. Delta file is input.

(firm ID) (delta)

* -f (str): final consumption output file (default null)

If the option is set, total final consumption of firms at every time step is output.

(time) (total final consumption)

* -F : fulldebug output

All internal variables are output at every time step. File name is fixed as "fullDebug.csv".

* -o (int): order algorithm type. 0: normal, 1: keep initial demand, 2: ignore negative inv adjustment

This option specifies which type of algorithm to be used to calculate order. Normal is the same as the original Hallegatte model. Keep initial demand is to ignore demand change and order the same level at the initial moment. Ignore negative inv(entory) adjustment is to ignore minus order and set the order to zero.

* -p (file:int:int...): get snapshot of production at indicated step(s). for more than one snapshot, int should be separated by ":"

You can get snapshot of PAct for all firms at designated time. This option requires the base filename and times when to take snapshots. For example, if

-p snapshot:10:20 

is set, you get snapshotAt10 and snapshotAt20 files.The format of file is

(firm ID) (PAct)

* -r (int): random seed (default 331)

To set random seed.

* -R (int): rationing algorithm. 0. proportional, 1. FC low priority, 2. lower has high priority, 3. 1 and 2 (default proportional)

This option can change rationing policy. It is indicated that the original rationing algorithm has a problem to simulate an actual network [Inoue and Todo, 2017].
	1. proportional: The original rationing policy. [Hallegatte, 2008; Henriet and Hallegatte, 2008]
	2. FC low priority: Final consumption is fulfilled last. (FC has low priority)
	3. lower has high priority: The improved rationing policy. [Inoue and Todo, 2017]
	4. 1 and 2: The algorithm is 2 but final consumption is fulfilled last in the same way as 1.

* -s (str): stats output file (default null)

This option specifies a stat filename to output simple stats of simulation. The following is output.

	1.totalPactutalStart
	2.totalValueAddedStart
	3.totalDamage
	4.totalDamageFirmNum
	5.totalPactutalEnd
	6.totalValueAddedEnd

* -t (int): simulation step (default 10)

To specify how many steps to be simulated.

* -v (str): value added output file (default null)

Total value added at each step is output.

(Time) (Total value added)

A value added is caluculated for each firm. The equation is (output)-(input). Minus value addded is acceptable.



To simulate the model in [Inoue and Todo, 2017], the following is the example of usage.

proNet.exe -R 2 -t 365 -v outputValueAddedToy.txt -d toyDelta01.txt toyATableLine.txt toyCVector.txt toyPiniVector.txt toyFirmAffiliation.txt outputPactToy.txt


Hope you can enjoy the code.
Comments are welcome.
