# Traveling Santa 2018 - Prime Paths (Kaggle Competition)#
This is the code I have created for the Traveling Santa 2018, my first Kaggle competition where I ranked 210/1874 (top 12%). This was a highly participated competition that featured top researchers like Bill Cook and Keld Helsgaun (who ended up winning the competition). The problem was very similar to a TSP with 197769 cities. In the figure below, the red node represents the first/last node.

![matlabmaptsp](https://user-images.githubusercontent.com/40466329/51650481-51ca5400-1f80-11e9-9c9f-fe4563d3aef1.jpg)

However, there was this twist:
The submission is scored on the Euclidean distance of your submitted path, subject to the constraint that every 10th step is 10% more lengthy unless coming from a prime CityId.
My approach was to first find the best Hamiltonian cycle disregarding the prime penalty and subsequently optimize from that initial solution.

## Step 1 - LKH ##
I used the [Lin-Kernighan-Helsgaun](http://akira.ruc.dk/~keld/research/LKH-3/) open-source code to find the best Hamiltonian cycle for this problem. After fine-tuning its parameters, I let it run for 2 days with occasional interruptions. In broad terms, my approach was:
1. To set a fairly long initial time to calculate the node penalties using gradient ascent in the first run (10 000 seconds)
2. On the first runs I used 8-opt sequential moves for the local search and progressively decreased this number until reaching 3-opt.
3. I used the maximum number of disjoint and alternating cycles to patch in order to find a gainful move. I tried relaxing the gain criterion (ie considering patches that had non-gainful edges) but the computation time increased substantially without a noticeable difference in the optimization.
4. Instead of a random walk, I used a double-bridge kick to avoid local minima, as most references I read mentioned this type of kick as yielding really good results without ruining the original path. 
5. Used the node candidates retrieved from Helsgaun's POPMUSIC metaheuristic, as numerous sources reported better results than using nearest neighbours or quadrant candidates. I defined the candidate's sample size to be 50 and the maximum number of neighbours to be 20 as per the POPMUSIC report of the LKH documentation 
6. Used a few initializations to avoid local minima caused by using the same seed for number generation and kept increasing the number of kicks for each initialization.
7. Multiplied the cities coordinates by 1e3 in order to get a higher precision in the distance calculation.

With this process I ended up with a pure distance solution that was 500 points better than the pure LKH solution shared in the kernels.

## Step 2 - Penalty optimization ##
The real challenge was this second part. Now that we have a fairly good tour, the idea is to create a gainful tradeoff between the increase of tour length to decrease the total penalty. For this, I implemented a 2.5opt and subsequently a 4+opt local search algorithm. For both algorithms, the main challenge was to code them in a way where they would not be too computationally expensive (as I was running this code on my laptop). For that, I used some code optimization strategies:
1.  Prime Set: Instead of constantly checking if a number is prime using MATLAB "isprime" function, I precomputed all the prime numbers and saved the prime number array. Moreover, to check if a CityID was in the prime array, I used the "ismembc" undocumented function that has a performance close to a C programming code.
2.  Permutations list: For both algorithms only a specific set of candidates is used for an opt move. The pool of candidates for each node is chosen based on the euclidean distance between a node and its neighbour. The search is done using a KDTree structure. All the possible node permutations are precomputed (see PermCalc.m)
3.  Edge cost precomputation for fast path scoring: When you perform a k-opt move, even when you just compute the cost of the path included in such move, there are many edges which cost remain unchanged from the original path. Capitalizing on this idea, I precomputed each edge distance from the original path and used this information to compute the cost of subpaths after a k-opt move. With this strategy, you do not have to compute the cost of the inter-node paths (even when they are reversed!) and only compute the new edges cost. The new penalties are also quite straightforward to compute (check code for details).

### 2.5Opt ###
This algorithm tries to find a gainful move by swapping 2 nodes. Moreover, it tries to reverse the path between the 2 nodes in order to find a gainful move. If the *maxPerms* variable is activated, the algorithm tests for moves across all the possible permutations between path and nodes like in the following example:

Be n1 and n2 the nodes to swap and p1 the path between them:

---- n1 - p1 - n2 ----

The 2.5-opt tests for all the 3! possible permutations of these 3 elements with b1 in normal and reversed order (total: 3 * 2^2 moves). File: TwoHalfOptV7.m

### 4+Opt ###
An extension of the previous algorithm with a few variables to control its computational cost. File: FourOptPlusV4.m

### Brute-force ###
When all the other methods were completed I ran a brute-force algorithm with O(n!) complexity for the final improvements. This usually resulted in improvements of 10-20 points. File: BruteForceFactorial.m

## Kernel Optimizations ##
The following Kernels represented, altogether, an improvement of about 60 points on my final solution and invaluable learning material as I am not yet very familiar with Python. Thanks to the authors!
* [Not a 5-and-5-halves-opt](https://www.kaggle.com/kostyaatarik/not-a-5-and-5-halves-opt) by *Kostya Atarik* which previous kernels were definitely an inspiration to create my own.
* [k-opt algorithm](https://www.kaggle.com/ubarredo/k-opt-algorithm) by Unai Barredo
* [DP Shuffle](https://www.kaggle.com/blacksix/dp-shuffle) by *blacksix* which was crucial to learn about Held-Karp's algorithm and dynamic programming.
* [Local Optimization using Google OR-Tool](https://www.kaggle.com/hblearn/local-optimization-using-google-or-tool) by *babachan* which was my first contact with Google's OR-tool.

## Bonus ##

### Simulated Annealing 3.5 ###
By the end of the competition, I tried a different kind of approach to local search algorithms with Simulated Annealing. In my implementation, the user can decide the number of initializations, runs ("reheats"), temperature and cooling schedule. I have also added a customizable tabu list to be able to induce more or less exploration in the annealing process. To the regular SA algorithm, I added a customizable 3.5-opt that runs when the temperature reaches 0. I was not able to fully test/optimize this implementation due to the large number of parameters to fine-tune. Therefore this code was not very useful for the competition. File: SimAnnOptV2.m

### Path Visualization ###
A script that allows a neat visualization of the map, tour and the penalty edges.

![tourvisual](https://user-images.githubusercontent.com/40466329/51650540-850ce300-1f80-11e9-903f-83518f237a23.jpg)

## Outstanding Ideas from high achievers in the competition ##
* Using GPX2 for tour recombination instead of IPT when doing the step 1;
* Using an LK or LKH method to deal with the penalty;
* Using POPMUSIC candidates for k-opt permutations in the step 2;
* Steadily incrementing the penalty instead of trying to optimize for the full 10% penalty right from the start.
* Doing more initializations/kicks.
