
package com.net2plan.general;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import com.jom.DoubleMatrixND;
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
		op.addConstraint("closest_n1n2(all,:) <= r_n"); // a closest node must be a node with a replica
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
	
	public static DoubleMatrix2D placeReplicasJavi (NetPlan netPlan, List<Node> dataCentersThisCDNThisApp , DoubleMatrix2D rtt_n1n2, int totalNumberOfReplicasToDistribute, double [] popularity, int U)
	{
		/* Basic checks */
//		final int N = netPlan.getNumberOfNodes();

		/* Create the optimization problem object (JOM library) */
		OptimizationProblem op = new OptimizationProblem();
		
		int numOfDCs = dataCentersThisCDNThisApp.size();
		int N = netPlan.getNumberOfNodes();
		
		/* Set some input parameters */
		op.setInputParameter("R", totalNumberOfReplicasToDistribute);
		
		System.out.println("num of DCs = " +  numOfDCs);		
		System.out.println("C = " + totalNumberOfReplicasToDistribute);
		
		double[][][] rtt_d = new double[U][N][numOfDCs];
		double[][][] p_u = new double[U][N][numOfDCs];
		
		for(int u = 0; u < U; u++)		
			for (int n = 0; n < N; n ++)			
				for (int d = 0; d < numOfDCs; d++) 	
				{
					p_u[u][n][d] = popularity[u];
					rtt_d[u][n][d] = rtt_n1n2.get(n, dataCentersThisCDNThisApp.get(d).getIndex());
				}
			
				
		op.setInputParameter("rtt_d", new DoubleMatrixND(rtt_d)); 												// mean RTT to the DC d from the rest of the network
		op.setInputParameter("p_u", new DoubleMatrixND(p_u));			// Popularity of the content units
		op.setInputParameter("CAPACITYDC", Math.ceil(totalNumberOfReplicasToDistribute / numOfDCs));			// Popularity of the content units
		
		/* Add the decision variables to the problem */
		op.addDecisionVariable("r_ud", true, new int[] { U, numOfDCs }, 0, 1); // there is a replica of u in dc d
		op.addDecisionVariable("closest_und", true, new int[] { U, N , numOfDCs }, 0, 1); // 1 if for cu u, the DC d is the closest to user in node n

		/* Sets the objective function */
		op.setObjectiveFunction("minimize", "sum (p_u .* closest_und .* rtt_d ) ");

		/* Constraints */
		for (int n = 0 ; n < N ; n ++)
		{
			op.setInputParameter ("n" , n);
			op.addConstraint("closest_und(all,n,all) <= r_ud "); // 
		}
		op.addConstraint("sum(r_ud) == R"); 					// we distribute all the replicas
		op.addConstraint("sum(r_ud,1) <= CAPACITYDC"); 			// no DC is oversubscribed
		op.addConstraint("sum(r_ud,2) >= 1"); 					// redundant

//		for (int n = 0 ; n < N ; n ++)
//		{
//			op.setInputParameter ("n" , n);
//			op.addConstraint("sum(closest_und(all,n,all)) == R "); 
//			op.addConstraint("sum(closest_und(all,n,all),2) >= 1");
//			op.addConstraint("sum(closest_und(all,n,all),1) <= CAPACITYDC");
//		}
//		for (int u = 0; u < U; u++)
//		{
//			op.setInputParameter ("u" , u);
//			op.addConstraint("sum(closest_und(u,all,all),2) == 1 ");
//		}
		
		
		/* Call the solver to solve the problem */
		op.solve("cplex", "solverLibraryName", "" );

		/* If no solution is found, quit */
		if (op.feasibleSolutionDoesNotExist()) throw new Net2PlanException("The problem has no feasible solution");
		if (!op.solutionIsFeasible()) throw new Net2PlanException("A feasible solution was not found");
		
		/* Retrieve the optimum solutions */
		DoubleMatrix2D  z_ud = op.getPrimalSolution("r_ud").view2D();
		for (int d1 = 0; d1 < numOfDCs; d1 ++)
			if (z_ud.viewColumn(d1).zSum() > totalNumberOfReplicasToDistribute) throw new RuntimeException();
		
		System.out.println("-------");
		for (int d1 = 0; d1 < rtt_d.length; d1++)
			System.out.println("Rtt("+d1+") = " + rtt_d[d1]);
		
		String res = "";
		for(int u = 0; u < U; u++)	
		{
			res += "U = " + u + " p_u = " + popularity[u] + " " ;
			for (int d2 = 0; d2 < numOfDCs; d2++)
				res += z_ud.get(u, d2) + " " ;
			System.out.println(res);
			res = "";
		}
		
		return z_ud;
	}
	
	
	
}