
package com.net2plan.general;

import java.util.List;

import com.jom.DoubleMatrixND;
import com.jom.OptimizationProblem;
import com.net2plan.interfaces.networkDesign.Net2PlanException;
import com.net2plan.interfaces.networkDesign.NetPlan;
import com.net2plan.interfaces.networkDesign.Node;
import com.net2plan.utils.Pair;

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
	
	
	public static Pair<DoubleMatrix2D,DoubleMatrixND> placeReplicas (NetPlan netPlan, List<Node> dataCentersThisCDNThisApp , DoubleMatrix2D rtt_n1n2, int totalNumberOfReplicasToDistribute, double [] popularity, double [] population_n , int U)
	{
		/* Basic checks */
//		final int N = netPlan.getNumberOfNodes();

		/* Create the optimization problem object (JOM library) */
		OptimizationProblem op = new OptimizationProblem();
		
		final int numOfDCs = dataCentersThisCDNThisApp.size();
		final int N = netPlan.getNumberOfNodes();
		final int capacityofDC = (int) Math.ceil(totalNumberOfReplicasToDistribute / (double) numOfDCs);
		/* Set some input parameters */
				
		double[][][] rtt_und = new double[U][N][numOfDCs];
		double[][][] popularity_und = new double[U][N][numOfDCs];
		double[][][] population_und = new double[U][N][numOfDCs];
		
		for(int u = 0; u < U; u++)		
			for (int n = 0; n < N; n ++)			
				for (int d = 0; d < numOfDCs; d++) 	
				{
					popularity_und[u][n][d] = popularity[u];
					population_und[u][n][d] = population_n[n];
					rtt_und[u][n][d] = rtt_n1n2.get(n, dataCentersThisCDNThisApp.get(d).getIndex());
				}			
				
		op.setInputParameter("R", totalNumberOfReplicasToDistribute);
		op.setInputParameter("rtt_und", new DoubleMatrixND(rtt_und)); 												// mean RTT to the DC d from the rest of the network
		op.setInputParameter("population_und", new DoubleMatrixND(population_und));			// Popularity of the content units
		op.setInputParameter("popularity_und", new DoubleMatrixND(popularity_und));			// Popularity of the content units
		op.setInputParameter("CAPACITYDC", capacityofDC);			// Popularity of the content units
		
		/* Add the decision variables to the problem */
		op.addDecisionVariable("r_ud", true, new int[] { U, numOfDCs }, 0, Integer.MAX_VALUE); // there is a replica of u in dc d
		op.addDecisionVariable("closest_und", true, new int[] { U, N , numOfDCs }, 0, 1); // 1 if for cu u, the DC d is the closest to user in node n

		/* Sets the objective function */
		op.setObjectiveFunction("minimize", "sum (popularity_und .* population_und .* closest_und .* rtt_und ) ");

		/* Constraints */
		op.addConstraint("sum(closest_und,3) == 1");
		for (int n = 0 ; n < N ; n ++)
		{
			op.setInputParameter ("n" , n);
			op.addConstraint("sum(closest_und(all,n,all),2) <= r_ud "); // 
		}
		op.addConstraint("sum(r_ud) == R"); 					// we distribute all the replicas
		op.addConstraint("sum(r_ud,1) <= CAPACITYDC"); 			// no DC is oversubscribed
		
		/* Call the solver to solve the problem */
		op.solve("cplex", "solverLibraryName", "" );

		/* If no solution is found, quit */
		if (op.feasibleSolutionDoesNotExist()) throw new Net2PlanException("The problem has no feasible solution");
		if (!op.solutionIsFeasible()) throw new Net2PlanException("A feasible solution was not found");
		
		/* Retrieve the optimum solutions */
		DoubleMatrix2D r_ud = op.getPrimalSolution("r_ud").view2D();
		DoubleMatrixND closest_und = op.getPrimalSolution("closest_und");
		
		for (int d = 0; d < numOfDCs; d ++)
			if (r_ud.viewColumn(d).zSum() > capacityofDC) throw new RuntimeException();
		for (int u = 0; u < U; u ++)
			if (r_ud.viewRow(u).zSum() == 0) throw new RuntimeException();
		if (r_ud.zSum() != totalNumberOfReplicasToDistribute) throw new RuntimeException();
		
//		String res = "";
//		for(int u = 0; u < U; u++)	
//		{
//			res += "";
//			for (int d = 0; d < numOfDCs; d++)
//				res += r_ud.get(u, d) + " " ;
//			System.out.println(res);
//			res = "";
//		}
//		for (int d = 0; d < numOfDCs; d++)
//			System.out.println("Num replicas of DC " + d + " = " + r_ud.viewColumn(d).zSum() );
		
		return new Pair<DoubleMatrix2D, DoubleMatrixND>(r_ud,closest_und, false);
	}
	
	
	
}