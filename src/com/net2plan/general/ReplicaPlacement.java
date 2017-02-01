package com.net2plan.general;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import com.jom.OptimizationProblem;
import com.net2plan.interfaces.networkDesign.Net2PlanException;
import com.net2plan.interfaces.networkDesign.NetPlan;
import com.net2plan.interfaces.networkDesign.Node;
import com.net2plan.utils.DoubleUtils;

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
		op.setInputParameter("c_n1n2", rtt_n1n2); // RTT resulting for the traffic of users in node n
		op.setInputParameter("K", numReplicas);
		final double [] u_n = new double [N]; 
		for (Entry<Node,Integer> entry : availableReplicaPositionsPerNode.entrySet()) 		
			u_n [entry.getKey().getIndex()] = entry.getValue();
		
		op.setInputParameter("u_n", u_n , "row");
		op.setInputParameter("pop_n", population_n , "row");
//		op.setInputParameter("M", 2E6);
		
		/* Add the decision variables to the problem */
		op.addDecisionVariable("r_n", true, new int[] { 1, N }, new double [N], u_n); // number of replicas in each node
		op.addDecisionVariable("closest_n1n2", true, new int[] { N, N }, 0, 1); 

		/* Sets the objective function */
		op.setObjectiveFunction("minimize", "sum (pop_n * (c_n1n2 .* closest_n1n2))");

		/* Constraints */
		op.addConstraint("sum(r_n) == K"); // distribute the given number of replicas
		op.addConstraint("sum(closest_n1n2,1) == 1"); // each node has exactly one closest one with a replica
//		op.addConstraint("closest_n1n2(all,:) <= r_n"); // a closest node must be a node with a replica
//		op.addConstraint("sum(closest_n1n2,2) <= r_n");
		
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
	
	public static DoubleMatrix2D placeReplicasJavi (NetPlan netPlan, List<Node> dataCentersThisCDNThisApp , DoubleMatrix2D rtt_n1n2, int capacityInNumberOfCU, double [] popularity, int U)
	{
		/* Basic checks */
//		final int N = netPlan.getNumberOfNodes();

		/* Create the optimization problem object (JOM library) */
		OptimizationProblem op = new OptimizationProblem();
		
		int numOfDCs = dataCentersThisCDNThisApp.size();

		
		/* Set some input parameters */
		op.setInputParameter("C", capacityInNumberOfCU);
		op.setInputParameter("D", numOfDCs);							// Number of DC this CDN
		op.setInputParameter("U", U);   								// Number of content units
		
		System.out.println("num of DCs = " +  numOfDCs);		
		System.out.println("C = " + capacityInNumberOfCU);
		
		double[] rtt_d = new double[numOfDCs];
		int d = 0;
		for (Node dc : dataCentersThisCDNThisApp) 	
		{
			rtt_d [d] = DoubleUtils.average(rtt_n1n2.viewColumn(dc.getIndex()).toArray());
			d++;
		}
		op.setInputParameter("rtt_d", rtt_d,"row"); 			// mean RTT to the DC d from the rest of the network
		op.setInputParameter("p_u", popularity, "row");			// Popularity of the content units
		
		/* Add the decision variables to the problem */
		op.addDecisionVariable("r_ud", true, new int[] { U, numOfDCs }, 0, 1); // number of replicas in each node

		/* Sets the objective function */
		op.setObjectiveFunction("maximize", "sum (p_u * r_ud * (1./rtt_d')) ");

		/* Constraints */
		op.addConstraint("sum(r_ud) <= 0.8*C"); // adjust the number of replicas to the 80% of the total capacity in the CDN 
		op.addConstraint("sum(r_ud,2) >= 1"); // At least one replica of each content unit in the CDN
		op.addConstraint("sum(r_ud,2) <= D"); // At most D replicas of each content unit in the CDN
		op.addConstraint("sum(r_ud,1) <= C"); // The content units must be at less than the capacity of the DC
		
		/* Call the solver to solve the problem */
		op.solve("cplex", "solverLibraryName", "" , "maxSolverTimeInSeconds" , 10.0);

		/* If no solution is found, quit */
		if (op.feasibleSolutionDoesNotExist()) throw new Net2PlanException("The problem has no feasible solution");
		if (!op.solutionIsFeasible()) throw new Net2PlanException("A feasible solution was not found");
		
		/* Retrieve the optimum solutions */
		DoubleMatrix2D  z_ud = op.getPrimalSolution("r_ud").view2D();
		for (int d1 = 0; d1 < numOfDCs; d1 ++)
			if (z_ud.viewColumn(d1).zSum() > capacityInNumberOfCU) throw new RuntimeException();
		
		System.out.println("-------");
		String res = "";
		for(int u = 0; u < U; u++)	
		{
			for (int d2 = 0; d2 < numOfDCs; d2++)
				res = res + z_ud.get(u, d2) + " " ;
			System.out.println(res);
			res = "";
		}
		
		return z_ud;
	}
	
	
	
}
