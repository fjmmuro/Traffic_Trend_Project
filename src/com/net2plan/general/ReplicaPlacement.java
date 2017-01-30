package com.net2plan.general;

import java.util.HashMap;
import java.util.Map;
import java.util.Map.Entry;

import com.jom.OptimizationProblem;
import com.net2plan.interfaces.networkDesign.Net2PlanException;
import com.net2plan.interfaces.networkDesign.NetPlan;
import com.net2plan.interfaces.networkDesign.Node;

import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;

/**
 * Solves several variants of node location problem formlations.
 * @net2plan.description
 * @net2plan.keywords Topology assignment (TA), JOM
 * @net2plan.ocnbooksections Section 7.2
 * @net2plan.inputParameters 
 * @author Pablo Pavon-Marino
 */
public class ReplicaPlacement 
{
	
	public static Map<Node,Integer> placeReplicas (NetPlan netPlan, Map<Node,Integer> availableReplicaPositionsPerNode , int numReplicas , DoubleMatrix2D rtt_n1n2 , double [] population_n)
	{
		/* Basic checks */
		final int N = netPlan.getNumberOfNodes();

		/* Create the optimization problem object (JOM library) */
		OptimizationProblem op = new OptimizationProblem();

		/* Set some input parameters */
		op.setInputParameter("c_n1n2", rtt_n1n2);
		op.setInputParameter("K", numReplicas);
		final double [] u_n = new double [N]; 
		for (Entry<Node,Integer> entry : availableReplicaPositionsPerNode.entrySet()) 
		{
			u_n [entry.getKey().getIndex()] = entry.getValue();
		}
		op.setInputParameter("u_n", u_n , "row");
		op.setInputParameter("pop_n", population_n , "row");
		op.setInputParameter("M", 2E6);
		
		/* Add the decision variables to the problem */
		op.addDecisionVariable("r_n", true, new int[] { 1, N }, new double [N], u_n); // number of replicas in each node
		op.addDecisionVariable("closest_n1n2", true, new int[] { N, N }, 0, 1); // RTT resulting for the traffic of users in node n

		/* Sets the objective function */
		op.setObjectiveFunction("minimize", "sum (pop_n * (c_n1n2 .* closest_n1n2))");

		/* Constraints */
		op.addConstraint("sum(r_n) == K"); // distribute the given number of replicas
		op.addConstraint("sum(closest_n1n2,1) == 1"); // each node has exactly one closest one with a replica
		op.addConstraint("closest_n1n2(all,:) <= r_n"); // a closest node must be a node with a replica
		
		/* Call the solver to solve the problem */
		op.solve("cplex", "solverLibraryName", "" , "maxSolverTimeInSeconds" , 10.0);

		/* If no solution is found, quit */
		if (op.feasibleSolutionDoesNotExist()) throw new Net2PlanException("The problem has no feasible solution");
		if (!op.solutionIsFeasible()) throw new Net2PlanException("A feasible solution was not found");
		
		/* Retrieve the optimum solutions */
		DoubleMatrix1D  z_j = op.getPrimalSolution("r_n").view1D();
		if (z_j.zSum() != numReplicas) throw new RuntimeException();
		Map<Node,Integer> numReplicasPerNode = new HashMap<> ();
		for (Node n : netPlan.getNodes())
		{
			final int numReplicasThisNode = (int) z_j.get(n.getIndex());
			if (numReplicasThisNode > availableReplicaPositionsPerNode.get(n)) throw new RuntimeException();
			numReplicasPerNode.put(n , numReplicasThisNode);
		}
		
		/* Check */
		
		return numReplicasPerNode;
	}
	
}
