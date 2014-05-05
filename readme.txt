= Tool for fitting Bayesian Knowledge Tracing models =

This tool is developed to fit Bayesian Knowledge Training models efficiently on
large datasets. It is a command line utility written in C/C++. Most of the
development was done while the author was a post-doc at Human-Computer
Interaction Institute, Carnegie Mellon University. These sources are published
under BSD-new (3-clause BSD) license.

= Bayesian Knowledge Tracing =

BKT is a user modeling approach in wide use in the area of Intelligent Tutoring
Systems. The goal of BKT is to infer whether the student has mastered a skill or
not from a pattern of successful or unsuccessful attempts to apply the skill. It
is a special case of Hidden Markov Model (HMM). In BKT there are two types of
nodes: binary state nodes capture skill mastery (mastered or not) and binary
observation nodes (correct or incorrect application of a skill).

For each skill BKT has four parameters.

1) pInit or pLo - is a probability the skill was known a priori,
2) pLearn or pT - is a probability the skill will transition into "mastered"
    state upon practice attempt,
3) pSlip or pS - is a probability that a known skill is applied incorrectly, and
4) pGuess or pG - is a probability that unknown skill will be applied correctly.

There is no forgetting in BKT and it's pForget or pF is set to zero. In
addition, there is a pKnown or pL parameter, which is a running estimate of the
skill mastery.

Parameters of the BKT can be represented in the matrix form. Priors -- \pi --
are a vector 1*N, where N is the number of states. Traditionally, for a 2-state
BKT, we assume that the first element of the \pi vector is the apriory
probability of knowing the skill. pLearn is part of the transitions matrix A
that captures probabilities of state chages from row to column and is N*N. No
forgetting is captured by setting A[1,2] = 0 (from mastered to unmastered).
pLearn corresponds to A[2,1] - from unmastered to mastered. pGuess and pSlip are
specified in B -- observation matrix N*M, where M is the number of observations.
First column corresponds to correct observation, second -- incorrect. For two
onservations, typical for BKT, pGuess is B[2,1] -- mastered skill but a correct
response, and pSlip is B[1,2] - mastered skill but incorrect response.

\pi .-------------------.
    |   pLo   | 1 - pLo |
    .-------------------.

A   .-------------------.
    |    1    |    0    |
    |-------------------|
    |    pT   | 1 - pT  |
    .-------------------.
    
B   .-------------------.
    | 1 - pS  |    pS   |
    |-------------------|
    |    pG   | 1 - pG  |
    .-------------------.

For more details on BKT refer to [1]. [2], among other things, discusses how
a gradient-based fitting of HMM can be implemented. [3, 4] covers additional
topics relevant for the implementation.

= Compilation =

If you are on a Linux or Mac OS X system, simply use 'make all' command. If you
are on Windows, you might need to install cygwin and have g++/gcc compiler
available and be sure to install 'make' command with cygwin.

= Data =

Input file data format is quite simple. Four tab separated columns: observation,
student, problem/problem step, skill(s). Observation is a 1-started integer. For
the two-state BKT, we advise to use 1 for 'correct' and 2 for 'incorrect'.
Student is a string label, so is problem or problem step, whatever granularity
you prefer. Skill is a string label. Multiple skill labels should be delimited
by a character of your choice (do not use tab). An example of few lines of input
is below where tilde symbol '~' is used as delimiter.

-- input file --
2   student_001 unit1-secion1-problem5-step1  addition~multiplication
1   student_001 unit1-secion1-problem5-step2  multiplication
1   student_001 unit1-secion1-problem5-step3  addition
-- input file --

If there is no skill label for a particular row of data use '.' (dot) symbol.
In test data, the utility will use known observations for training and will
produce predictions for missing observations that should have '.' (dot) instead
of observation code (1, 2, or otherwise).

Output files are the model and the predictions file. Model file contains
general information about the data, e.g., number of observations, number of
skills and number of problems/problem steps; as well as the model parameters.

Prediction file consists of model predictions for each row of the input file.
Depending on the option (see parameter specifications below), you can print out
probability distribution for student's response (probability of correct and
probability of incorrect, if the observation node is binarye) and, additionally,
print out probability distributions over the values of the hidden state for all
skills specified for the data point. The probability values are tab separated.
See examples of a outputs file below.

-- model file --
SolverId	1.1
nK	1
nG	1
nS	2
nO	2
Null skill ratios	  1.0000000	  0.0000000
0	multiplication-skill
PI	0.50000000	0.50000000
A	1.00000000  0.00000000	0.40000000	0.60000000
B	0.80000000	0.20000000	0.20000000	0.80000000
-- model file --

In the model file example above, we see a specification of the solver algorithm
(to be discussed below), one skill (nK=1), one student in the data (nG=1), two
hidden states and 2 observations (nS and nO respectively). Skill numbering is
zero started. Skill name is "multiplication-skill". Rows prefixed with PI, A,
and B correspond to priors, transition, and observation matrices written row by
row. Thus, first element of PI is pLo, 3rd element of A is pT, 2nd and 3rd
elements of B are pS and pG respectively.

-- prediction file --
0.73    0.27
0.88    0.12
0.94    0.06
0.99    0.01
-- prediction file --


= Training BKT models =

trainhmm utility has the following launch signature:

trainhmm [options] input_file [[model_file] prediction_file]
options:
-s : structure.solver[.solver setting], structures: 1-by skill, 2-by user;
     solvers: 1-Baum-Welch, 2-Gradient Descent, 3-Conjugate Gradient Descent;
     Conjugate Gradient Descent has 3 settings: 1-Polak-Ribiere,
     2-Fletcher–Reeves, 3-Hestenes-Stiefel.
     For example '-s 1.3.1' would be by skill structure (classical) with
     Conjugate Gradient Descent and Hestenes-Stiefel formula, '-s 2.1' would be
     by student structure fit using Baum-Welch method.
-S : perform scaling of forward/backward variables: 0 - off (default), 1 - on.
     Only allowed for Baum-Welch solver ('-s 1.1' setting), otherwise auto-set
     to off.
-e : tolerance of termination criterion (0.01 default)
-i : maximum number of iterations (200 by default)
-q : quiet mode, without output, 0-no (default), or 1-yes
-n : number of hidden states, should be 2 or more (default 2)
-0 : initial parameters comma-separated for priors, transition, and emission
     probabilities skipping the last value from each vector (matrix row) since
     they should sum up to 1; default 0.5,1.0,0.4,0.8,0.2
-l : lower boundaries for parameters, comma-separated for priors, transition,
     and emission probabilities (without skips); default 0,0,1,0,0,0,0,0,0,0
-u : upper boundaries for params, comma-separated for priors, transition,
     and emission probabilities (without skips); default 0,0,1,0,0,0,0,0,0,0
-c : weight of the L2 penalty, 0 (default)
-f : fit as one skill, 0-no (default), 1 - fit as one skill and use params as
     starting point for multi-skill, 2 - force one skill
-m : report model fitting metrics (AIC, BIC, RMSE) 0-no (default), 1-yes. To 
     specify observation for which metrics to be reported, list it after ','.
     For example '-m 0', '-m 1' (by default, observation 1 is assumed), '-m 1,2'
     (compute metrics for observation 2). Incompatible with '-v' option.
-v : cross-validation folds, stratification, and target state to validate
     against, default 0 (no cross-validation),
     examples '-v 5,i,2' - 5 fold, item-stratified c.-v., predict state 2,
     '-v 10' - 10-fold subject-stratified c.-v. predict state 1 by default,
     alternatively '-v 10,g,1', and finally '-v 5,n,2,' - 5-fold unstratified
     c.-v. predicting state 1.
-p : report model predictions on the train set 0-no (default), 1-yes; 2-yes,
     plus output state probability; works with -v and -m parameters.
-d : delimiter for multiple skills per observation; 0-single skill per
     observation (default), otherwise -- delimiter character, e.g. '-d ~'.
-b : treat input file as binary input file created from text file by
     inputconvert utility.
-B : block re-estimation of prior, transitions, or emissions parameters
     respectively (defailt is '-B 0,0,0'), to block re-estimation of transition
     probabilities specify '-B 0,1,0'.


= Using models for prediction =

predicthmm utility has the following launch signature:

predicthmm [options] input_file model_file [predicted_response_file]
options:
-q : quiet mode, without output, 0-no (default), or 1-yes
-m : report model fitting metrics (AIC, BIC, RMSE) 0-no (default), 1-yes. To 
     specify observation for which metrics to be reported, list it after ','.
     For example '-m 0', '-m 1' (by default, observation 1 is assumed), '-m 1,2'
     (compute metrics for observation 2). Incompatible with-v option.
-d : delimiter for multiple skills per observation; 0-single skill per
     observation (default), otherwise -- delimiter character, e.g. '-d ~'.
-b : treat input file as binary input file  created from text file by
     inputconvert utility.

= Examples =

Small sample data file <toy_data.txt> is generated using the following BKT
parameters: pLo=0.4, pT=0.35, pS=0.25, pG=0.12 .

To fit a BKT model of this data using an EM algorithm run the following command:

sh> ./trainhmm -s 1.1 -m 1 -p 1 toy_data.txt model.txt predict.txt

The model will have 90% accuracy and root mean squared error (RMSE) = 0.3227 and
the recovered BKT parameters would be: pLo=0.50, pT=0.17, pS=0.00, pG=0.00 .

If we fit BKT model using Conjugate Gradient method using '-s 1.2' argument, the
recovered parameters would be: pLo=0.00, pT=0.18, pS=0.08, pG=0.03, the accuracy
would remain at 90% while RMSE = 0.2982.

To generate predictions using a previously fit model run the following command 
(do not forget that prediction will only be generated for rows where observation
is now known -- marked with '.'):

sh> ./predicthmm toy_data_test.txt model.txt predict.txt

To give this tool a proper test you might want to try it on a KDD Cup 2010
dataset donated to the Pittsburgh Science of Learning Center by Carnegie
Learning Inc. The dataset can be downloaded (after a quick registration) from
http://pslcdatashop.web.cmu.edu/KDDCup/. This datasets consists of training and
challenge sets. For the sake of testing the tool, download the challenge 
Algebra I set that has about 9 million transactions of over 3300 students. The
training file should be trimmed to the tool's format. See shell commands below
that do that.

sh> gawk '-F\t' 'BEGIN{OFS=""} {print ".","\t",$2,"\t",$3,"__",$4,"\t",tolower($20)}' algebra_2008_2009_train.txt > tmp1.txt
sh> sed 1d tmp1.txt
sh> rm tmp1.txt
sh> awk '-F\t' 'BEGIN{OFS=""} {print $1,"\t",$2,"\t",$3,"\t",((length($4)==0)?".":$4)}' tmp2.txt > a89_kts_train.txt
sh> rm tmp2.txt

To fit a BKT model of this dataset using gradient descent method as well as to 
compute fit metrics and the prediction run the following command:

sh> ./trainhmm -s 1.2 -m 1 -p 1 a89_kts_train.txt model.txt predict.txt

Depending on your hardware, the model should be fit in approximately 2 minutes.

= References =

[1] Corbett, A. T. and Anderson, J. R.: Knowledge tracing: Modeling the
    acquisition of procedural knowledge. User Modeling and User-Adapted
    Interaction, 4(4), 253-278. (1995)
[2] Levinson, S. E., Rabiner, L. R., and Sondhi, M. M.: An Introduction to the
    Application of the Theory of Probabilistic Functions of a Markov Process to
    Automatic Speech Recognition. Bell System Technical Journal, 62(4):
    1035-1074. (1983)
[3] http://en.wikipedia.org/wiki/Wolfe_conditions
[4] http://en.wikipedia.org/wiki/Nonlinear_conjugate_gradient_method

= Contact Us =

If you have any questions, comments, propositions, feel free to contact me at
myudelson@gmail.com

Have fun fitting BKT models,
Michael (Mikhail) Yudelson
