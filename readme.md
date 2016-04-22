# Tool for fitting Bayesian Knowledge Tracing models

This tool is developed to fit Bayesian Knowledge Training models on large datasets command-line utility for fitting standard BKT models at scale. It is a command line utility written in C/C++. The project started while the author was a post-doc at Human-Computer Interaction Institute, Carnegie Mellon University. These sources are published under BSD-new (3-clause BSD) license.

# Bayesian Knowledge Tracing

BKT is a user modeling approach widely use in the area of Intelligent Tutoring Systems (ITS). In ITS, it is customary to tag problems and problem steps students are solving with knowledge quanta often called skills. When a student is attempting to solve a problem or a problem step, they are said to _practice a skill_. A student model in a general sense performs deductions on whether the skill(s) is (are) mastered given a sequence of correct and incorrect attempts to solve a problems/problem steps that are tagged with it. The goal of BKT is to infer whether the student has mastered a skill or not from a pattern of successful or unsuccessful attempts to apply it. BKT uses a formalism of a Hidden Markov Model to estimate student's skill mastery. In BKT, there are two types of nodes: binary state nodes capture skill mastery (one per skill) that assume values of _mastered_ or _not mastered_ and binary observation nodes (one per skill) that assume values of _correct_ or _incorrect_. 

Each skill in standard BKT has four parameters.

1. _p-init_ or p(Lo) - is a probability the skill was known a priori,
2. _p-learn_ or p(T) - is a probability the skill will transition into _mastered_ state upon practice attempt,
3. _p-forget_ or p(F) - is a probability that the skill will transition into _not mastered state_ upon practice attempt. Traditionally, _p-forget_ is set to zero and is not counted towards the total number of parameters.
4. _p-slip_ or p(S) - is a probability that a known skill is applied incorrectly, and
5. _p-guess_ or p(G) - is a probability that unknown skill will be applied correctly.

Parameters of the BKT can be represented in the matrix form. Priors – ![\pi](http://latex.codecogs.com/gif.latex?%5CPi) – are a vector 1\*N, where N is the number of states (N=2 in standard BKT). Throughout this discussion, we will we assume that the first element of the ![\pi](http://latex.codecogs.com/gif.latex?%5CPi) vector is the a priory probability of knowing the skill. _p-learn_ is part of the transitions matrix A that captures probabilities of state changes from the row to the column and is N\*N. No forgetting is captured by setting A[1,2] = 0 (from mastered to unmastered). _p-learn_ corresponds to A[2,1] - from unmastered to mastered. _p-guess_ and _p-slip_ are specified in B – observation matrix N\*M, where M is the number of observations (M=2 in standard BKT). The first column corresponds to a correct observation, the second – to the incorrect. For standard BKT with two observations, _p-guess_ is B[2,1] – unmastered skill but a correct response, and _p-splip_ is B[1,2] - mastered skill but an incorrect response. 

    
![\Pi=\[p(L_0), 1-p(L_0)\]](http://latex.codecogs.com/gif.latex?%5CPi%3D%0D%0A%5Cbegin%7Bbmatrix%7D%0D%0Ap%28L_0%29%261-p%28L_0%29%5C%5C%0D%0A%5Cend%7Bbmatrix%7D)

![A=\[1, 0; p(T), 1-p(T)\]](http://latex.codecogs.com/gif.latex?A%3D%0D%0A%5Cbegin%7Bbmatrix%7D%0D%0A1%260%5C%5C%0D%0Ap%28T%29%261-p%28T%29%0D%0A%5Cend%7Bbmatrix%7D)

![B=\[1-p(S), p(S); p(G), 1-p(G)\]](http://latex.codecogs.com/gif.latex?B%3D%0D%0A%5Cbegin%7Bbmatrix%7D%0D%0A1-p%28S%29%26p%28S%29%5C%5C%0D%0Ap%28G%29%261-p%28G%29%0D%0A%5Cend%7Bbmatrix%7D)


For more details on BKT refer to [1]. [2], among other things, discusses how a gradient-based fitting of HMM can be implemented. [3, 4] cover additional topics relevant for the implementation.

# Getting and Compiling the Code

A public version of the standard BKT code is available via International Educational Data Mining Society GitHub [repository](http://goo.gl/5DqpkW). If you have `git` utility installed, use the following command to get the source code.

```sh
git clone https://github.com/IEDMS/standard-bkt
```

To compile issue `make all` command. If you are on Linux, the `g++/gcc` compiler and Open MP library should already be installed. On Mac OX, you would need command line tools of Xcode to be installed. `g++/gcc` compiler with Open MP (no longer bundled with Mac OS X by default) could be downloaded from [hpc.sourceforge.net](http://hpc.sourceforge.net). If you are on Windows, you might need to install `cygwin` and have `g++/gcc` compiler available and be sure to install `make` command with `cygwin`.

# Data 

Input file data format is quite simple – four tab-separated columns: observation, student, problem/problem step, skill(s). Observation is a 1-started integer. For the two-state BKT, we advise using 1 for _correct_ and 2 for _incorrect_. The student is a string label, so is problem or problem step, whatever granularity you prefer. Skill is a string label. Multiple skill labels should be delimited by a character of your choice (do not use tab). An example with a few lines of input is below. Here, tilde symbol `~` is used as a delimiter. 

```
-- input file --
2   student_001 unit1-section1-problem5-step1  addition~multiplication
1   student_001 unit1-section1-problem5-step2  multiplication
1   student_001 unit1-section1-problem5-step3  addition
-- input file --
```

If there is no skill label for a particular row of data use `.` (dot) symbol. In test data, the utility will use known observations for training and will produce predictions for missing observations that should have `.` (dot) instead of observation code (1, 2, or otherwise).

Output files are the model file and the predictions file. The model file contains general information about the data, e.g., the number of observations, the number of skills and the number of problems/problem steps; as well as the model parameters.

Prediction file consists of model predictions for each row of the input file. Depending on the option (see parameter specifications below). You can print out probability distribution for student's response (probability of correct and probability of incorrect, if the observation node is binary) and, additionally, print out probability distributions over the values of the hidden state for all skills specified for that particular data point. The probability values are tab separated. See examples of a outputs file below. 

```
-- model file --
SolverId	1.1
nK	1
nG	1
nS	2
nO	2
Null skill ratios	  1.0000000	  0.0000000
0	multiplication-skill
PI	0.45000000	0.55000000
A	1.00000000  0.00000000	0.40000000	0.60000000
B	0.79000000	0.21000000	0.22000000	0.78000000
-- model file --
```

In the model file example above, we see a specification of the solver algorithm (to be discussed below), one skill (nK=1), one student in the data (nG=1), two hidden states and 2 observations (nS and nO respectively). Skill numbering is zero-started. Skill name is "multiplication-skill". Rows prefixed with PI, A, and B correspond to priors, transition, and observation matrices written row by row. Thus, the first element of PI is p(Lo)=0.55, 3rd element of A is p(T)=0.4, 2nd and 3rd elements of B are p(S)=0.12 and p(G)=0.22 respectively.

```
-- prediction file --
0.73    0.27
0.88    0.12
0.94    0.06
0.99    0.01
-- prediction file --
```

# Training BKT models 

```trainhmm``` utility has the following launch signature:

```
trainhmm [options] input_file [[model_file] prediction_file]
```

options:

`-s` : structure.solver[.solver setting], structures: 1-by skill, 2-by user; solvers: 1-Baum-Welch, 2-Gradient Descent, 3-Conjugate Gradient Descent; Conjugate Gradient Descent has 3 settings: 1-Polak-Ribiere, 2-Fletcher–Reeves, 3-Hestenes-Stiefel, and 4-Dai-Yuan. For example `-s 1.3.1` would be by skill structure (classical) with Conjugate Gradient Descent and Hestenes-Stiefel formula, `-s 2.1` would be by student structure fit using Baum-Welch method.

`-S` : perform scaling of forward/backward variables: 0 - off (default), 1 - on. Only allowed for Baum-Welch solver (`-s 1.1` setting), otherwise auto-set to off.

`-e` : tolerance of termination criterion (0.01 for parameter change default); could be computed by the change in log-likelihood per datapoint, e.g. `-e 0.00001,l`.

`-i` : maximum number of iterations (200 by default)

`-q` : quiet mode, without output, 0-no (default), or 1-yes

`-n` : number of hidden states, should be 2 or more (default is 2)

`-0` : initial parameters comma-separated for priors, transition, and emission probabilities skipping the last value from each vector (matrix row) since they should sum up to 1; default is `0.5,1.0,0.4,0.8,0.2`.

`-l` : lower boundaries for parameters, comma-separated for priors, transition, and emission probabilities (without skips); default is `0,0,1,0,0,0,0,0,0,0`.

`-u` : upper boundaries for params, comma-separated for priors, transition, and emission probabilities (without skips); default is `1,1,1,1,1,1,1,0.3,0.3,1` with slip and guess capped at 0.3.

`-c` : specification of the C weight and centroids for L2 penalty, empty (default). For standard BKT - 4 comma-separated numbers: C weight of the penalty and centroids, for PI, A, and B matrices respectively. If used for iBKT with student effects, 8 values will be used with 4 additional values for student effect matrices. For example, `-c 1.0,0.5,0.5,0.0`.

`-f` : fit as one skill, 0-no (default), 1 - fit as one skill and use params as starting point for multi-skill, 2 - force one skill

`-m` : report model fitting metrics (AIC, BIC, RMSE) 0-no (default), 1-yes. To  specify observation for which metrics to be reported, list it after `,`. For example `-m 0`, `-m 1` (by default, observation 1 is assumed), `-m 1,2` (compute metrics for observation 2). Incompatible with `-v` option.

`-v` : cross-validation folds, stratification, and target state to validate against, folds input/output file, default 0 (no cross-validation), examples `-v 5,i,2` - 5 fold, item-stratified c.-v., predict state 2, `-v 10` - 10-fold subject-stratified c.-v. predict state 1 by default, alternatively `-v 10,g,1`, `-v 5,n,2,folds.txt,o` - 5-fold unstratified c.-v. predicting state 2, [o]output folds to `folds.txt`, and here  `-v 5,n,2,folds.txt,i`, folds are actually read [i]n from the file.

`-p` : report model predictions on the train set 0-no (default), 1-yes; 2-yes, plus output state probability; works with -v and -m parameters.

`-U` : controls how the probability distribution of the states is updated. Takes the following format `-U r|g[,t|g]`, where the first character controls how prediction treats known observations, second – how prediction treats unknown observations, and third – whether to output probabilities of priors. Dealing with known observations `r` - reveal actual observations for the update of state probability distribution (makes sense for modeling how an actual system would work), `g` – `guessing` the observation based on the predicted outcomes (arg max) – more appropriate when comparing models (so that no information about observation is never revealed). Dealing with unknown observations (marked as `.` – dot): `t` – use transition matrix only, `g` – `guess` the observation. Default (if ommitted) is `-U r,t`. For examle, `-U g,g` would require _guessing_ of what the observation was using model parameters and the running value of the probabilities of state distributions.

`-d` : delimiter for multiple skills per observation; single skill per observation (default), otherwise – delimiter character, e.g. `-d ~`.

`-b` : treat input file as binary input file created from the text file by the `inputconvert` utility.

`-B` : block re-estimation of prior, transitions, or emissions parameters respectively (default is `-B 0,0,0`), to block re-estimation of transition probabilities specify `-B 0,1,0`.

`-P` : use parallel processing, default - 0 (no parallel processing), 1 - fit separate skills/students separately, 2 - fit separate sequences within skill/student separately.

`-o` : in addition to printing to the console, print output to the file specified after this option (by default is empty and no extra printing is done).


# Using models for prediction

```predicthmm``` utility has the following launch signature:

```
predicthmm [options] input_file model_file [predicted_response_file]
```

options:

`-q` : quiet mode, without output, 0-no (default), or 1-yes

`-m` : report model fitting metrics (AIC, BIC, RMSE) 0-no (default), 1-yes. To  specify observation for which metrics to be reported, list it after `,`. For example `-m 0`, `-m 1` (by default, observation 1 is assumed), `-m 1,2` (compute metrics for observation 2). Incompatible with-v option.

`-d` : delimiter for multiple skills per observation; single skill per observation (default), otherwise – delimiter character, e.g. `-d ~`.

`-b` : treat input file as binary input file created from a text file by the `inputconvert` utility.

`-p` : report model predictions on the train set 0-no (default), 1-yes; 2-yes, plus output state probability; works with `-v` and `-m` parameters.

`-U` : controls how the probability distribution of the states is updated. Takes the following format `-U r|g[,t|g]`, where the first character controls how prediction treats known observations, second – how prediction treats unknown observations, and third – whether to output probabilities of priors. Dealing with known observations `r` - reveal actual observations for the update of state probability distribution (makes sense for modeling how an actual system would work), `g` - `guessing` the observation based on the predicted outcomes (arg max) – more appropriate when comparing models (so that no information about observation is never revealed). Dealing with unknown observations (marked as `.` – dot): `t` – use transition matrix only, `g` – `guess` the observation. Default (if ommitted) is `-U r,t`. For examle, `-U g,g` would require `guessing` of what the observation was using model parameters and the running value of the probabilities of state distributions.


# Examples

Small sample data file <toy_data.txt> is generated using the following BKT parameters: pLo=0.4, pT=0.35, pS=0.25, pG=0.12. 

To fit a BKT model of this data using an EM algorithm run the following command:

```sh
./trainhmm -s 1.1 -m 1 -p 1 toy_data.txt model.txt predict.txt
```

The model will have 90% accuracy and root mean squared error (RMSE) = 0.302691 and the recovered BKT parameters would be: pLo=0.00000000, pT=0.16676161, pS=0.00044059, pG=0.00038573. Overall loglikelihood, actually, goes up from 9.3763477 to 10.4379501 in 3 iterations.

If we fit BKT model using Gradient Descent method using `-s 1.2` argument, the recovered parameters would be: pLo=0.00041944, pT=0.17478539, pS=0.07938036, 0.03804388, the accuracy would remain at 90% while RMSE = 0.299250. Loglikelihood changes from  9.3763477 to 6.4099682 after 11 iterations.

To generate predictions using a previously fit model run the following command  (do not forget that prediction will only be generated for rows where observation is not known – marked with `.`): 

```sh
./predicthmm -p 1 toy_data_test.txt model.txt predict.txt
```

To give this tool a proper test you might want to try it on a KDD Cup 2010 dataset donated to the Pittsburgh Science of Learning Center by Carnegie Learning Inc. The dataset can be downloaded (after a quick registration) from [here](http://pslcdatashop.web.cmu.edu/KDDCup/). This dataset consists of training and challenge sets. For the sake of testing the tool, download the challenge  Algebra I set that has about 9 million transactions of over 3300 students. The training file should be trimmed to the tool's format. See shell commands below that do that.

```sh
gawk -F"\t" 'BEGIN{OFS=""} {if(NR==1)next; skill=$20; gsub("~~", "~", skill); skill=(skill=="")?".":($3"__"skill); print 2-$14,$2,$3"__"$4,skill;}' algebra_2008_2009_train.txt > a89_kts_train.txt
```

To fit a BKT model of this dataset using gradient descent method as well as to compute fit metrics and the prediction run the following command:

```sh
./trainhmm -s 1.2 -d ~ -m 1 -p 1 a89_kts_train.txt model.txt predict.txt
```

Depending on your hardware, the model should be fit in about 1-2 minutes.

# References

1. Corbett, A. T., and Anderson, J. R.: Knowledge tracing: Modeling the acquisition of procedural knowledge. User Modeling and User-Adapted Interaction, 4(4), 253-278. (1995)
 
2. Levinson, S. E., Rabiner, L. R., and Sondhi, M. M.: An Introduction to the Application of the Theory of Probabilistic Functions of a Markov Process to Automatic Speech Recognition. Bell System Technical Journal, 62(4): 1035-1074. (1983)

3. Wolfe conditions. (2016, January 20). In _Wikipedia, The Free Encyclopedia_. [http://en.wikipedia.org/wiki/Wolfe_conditions](http://en.wikipedia.org/wiki/Wolfe_conditions)

4. Nonlinear conjugate gradient method. (2016, February 22). In _Wikipedia, The Free Encyclopedia_. [http://en.wikipedia.org/wiki/Nonlinear_conjugate_gradient_method](http://en.wikipedia.org/wiki/Nonlinear_conjugate_gradient_method)

# Contact Us

If you have any questions, comments, propositions, feel free to contact me at myudelson@gmail.com.

Have fun fitting BKT models,

Michael (Mikhail) Yudelson