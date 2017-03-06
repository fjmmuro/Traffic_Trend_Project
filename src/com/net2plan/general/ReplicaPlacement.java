package com.net2plan.general;

import java.util.*;

import cern.colt.matrix.tdouble.impl.DenseDoubleMatrix1D;
import com.google.common.collect.Sets;
import com.jom.DoubleMatrixND;
import com.jom.OptimizationProblem;
import com.net2plan.interfaces.networkDesign.Net2PlanException;
import com.net2plan.interfaces.networkDesign.NetPlan;
import com.net2plan.interfaces.networkDesign.Node;
import com.net2plan.utils.DoubleUtils;
import com.net2plan.utils.Pair;
import com.net2plan.utils.Triple;

import cern.colt.matrix.tdouble.DoubleFactory1D;
import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;

import javax.sound.midi.SysexMessage;

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
	
	/**
	 * @param netPlan
	 * @param previousDCPositionsThisCDN
	 * @param cost_n1n2
	 * @param popularity_u
	 * @param population_n
	 * @param h_a weight of application a in total traffic
	 * @param beta_a ration between user-DC vs Dc-DC traffic of application a
	 * @param path
	 * @return
	 */
	public static Pair<DoubleMatrix2D,List<Node>> placeReplicas (NetPlan netPlan, List<Node> previousDCPositionsThisCDN , 
			DoubleMatrix2D cost_n1n2, double [] popularity_u, double [] population_n , 
			DoubleMatrix1D h_a , DoubleMatrix1D beta_a , 
			String path, boolean isFirstTime)
	{
		Triple<DoubleMatrix2D , List<Node> , Double> bestSolution = null; // r_ud, dcPositions , objectiveCost
		final Set<Node> potentialNewDcPositions = Sets.difference(new HashSet<> (netPlan.getNodes()), new HashSet<> (previousDCPositionsThisCDN));
		double currentSolutionCost = Double.MIN_VALUE;

		boolean flag = true;
		do {
			for (Node newDCNode : potentialNewDcPositions)
			{
				List<Node> dcPositions = new ArrayList<>(previousDCPositionsThisCDN);
				if (!isFirstTime) dcPositions.add(newDCNode);

				/* Create the optimization problem object (JOM library) */
				OptimizationProblem op = new OptimizationProblem();
				final int U = popularity_u.length;
				final int numOfDCs = dcPositions.size();
				final int N = netPlan.getNumberOfNodes();
				final int A = (int) h_a.size();

	//			final int capacityofDC = (int) Math.ceil(totalNumberOfReplicasToDistribute / (double) oldNumOfDCs);
				int solverTime = 60;
				if (N > 20)
					solverTime = 150;

//				System.out.println("----- Popularity-----------");
//				System.out.println(Arrays.toString(popularity_u));

				/* Set some input parameters */
				double[][][] cost_uand = new double[U * A][N][numOfDCs];
				double[][][] trafUserToDC_uand = new double[U * A][N][numOfDCs];
				double[][] trafUserToDC_uad = new double[U * A][numOfDCs];
				double[][] beta_uad = new double[U * A][numOfDCs];
//				double[] e_d = new double[N];

				for (int u = 0; u < U; u++)
					for (int a = 0; a < A; a++)
						for (int n = 0; n < N; n++)
							for (int d = 0; d < numOfDCs; d++)
							{
								trafUserToDC_uand[u + a * U][n][d] = popularity_u[u] * population_n[n] * h_a.get(a) / DoubleUtils.sum(population_n);
//								beta_uad[u + a * U][d] = beta_a.get(a);
								cost_uand[u + a * U][n][d] = cost_n1n2.get(n, dcPositions.get(d).getIndex());
								trafUserToDC_uad[u + a * U][d] = popularity_u[u] * h_a.get(a) * beta_a.get(a);
//                                trafUserToDC_uad[u + a * U][d] = popularity_u[u] * population_n[n] * h_a.get(a);
							}

				System.out.println("A : " + A);

	//			op.setInputParameter("R", totalNumberOfReplicasToDistribute);
				op.setInputParameter("cost_uand", new DoubleMatrixND(cost_uand));                                                // mean RTT to the DC d from the rest of the network
				op.setInputParameter("trafUserToDC_uand", new DoubleMatrixND(trafUserToDC_uand));            // Popularity of the content units
				op.setInputParameter("trafUserToDC_uad", new DoubleMatrixND(trafUserToDC_uad));            // Popularity of the content units
//				op.setInputParameter("beta_uad", new DoubleMatrixND(beta_uad));           				 // Popularity of the content units
                op.setInputParameter("R",U*A);
				//			op.setInputParameter("numOfDCs",numOfDCs);
				//			op.setInputParameter("CAPACITYDC", capacityofDC);			// Popularity of the content units
				//			DoubleMatrix1D alreadyHasADc = DoubleFactory1D.dense.make(N , 0.0); for (Node n : previousDCPositionsThisCDN) alreadyHasADc.set(n.getIndex() , 1);

				/* Add the decision variables to the problem */
				op.addDecisionVariable("r_uad", true, new int[]{U * A, numOfDCs}, 0, 1); // there is a replica of u in dc d
				op.addDecisionVariable("closest_uand", true, new int[]{U * A, N, numOfDCs}, 0, 1); // 1 if for cu u, the DC d is the closest to user in node n

				/* Sets the objective function */
//				op.setObjectiveFunction("minimize", "sum (trafUserToDC_uand .* closest_uand .* cost_uand )"
//						+ " +   sum (trafUserToDC_uad' * (sum(r_uad,2) - 1)) ");

                op.setObjectiveFunction("minimize", "sum (trafUserToDC_uand .* closest_uand .* cost_uand )");
//                op.setObjectiveFunction("minimize", "sum (trafUserToDC_uad' * (sum(r_uad,2) - 1))");

				/* Constraints */
				op.addConstraint("sum(closest_uand,3) == 1");
				op.addConstraint("sum(r_uad,2) >= 2");
				op.addConstraint("sum(r_uad,1) >= 0.1*R");

				for (int n = 0; n < N; n++) {
					op.setInputParameter("n", n);
					op.addConstraint("sum(closest_uand(all,n,all),2) <= r_uad "); //
				}
							op.addConstraint("sum(r_uad) == 2*R"); 					// we distribute all the replicas

				/* Call the solver to solve the problem */
				op.solve("cplex", "solverLibraryName", path, "maxSolverTimeInSeconds", solverTime);

				/* If no solution is found, quit */
				if (op.feasibleSolutionDoesNotExist())
					throw new Net2PlanException("The problem has no feasible solution");
				if (!op.solutionIsFeasible()) throw new Net2PlanException("A feasible solution was not found");

				/* Retrieve the optimum solutions */
				DoubleMatrix2D r_uad = op.getPrimalSolution("r_uad").view2D();
				System.out.println(" Optimal cost: "+ op.getOptimalCost());
//				System.out.println("-------------------");
//				System.out.println(r_uad);


				System.out.println(r_uad.zMult(DoubleFactory1D.dense.make(numOfDCs,1.0),null));

				System.out.println(r_uad.viewDice().zMult(DoubleFactory1D.dense.make(U*A,1.0),null));


				//			for (int d = 0; d < oldNumOfDCs; d ++)
				//				if (r_ud.viewColumn(d).zSum() > capacityofDC) throw new RuntimeException();
				for (int u = 0; u < U; u++)
					if (r_uad.viewRow(u).zSum() < 2) throw new RuntimeException("At least one content unit has not two replicas");
				//			if (r_ud.zSum() != totalNumberOfReplicasToDistribute) throw new RuntimeException();
				for (int d = 0; d < numOfDCs; d++)
					if ((r_uad.viewColumn(d).zSum() == 0)) throw new RuntimeException("At least one DC is empty");
				//			if (e_d.zSum() != numOfDCs) throw new RuntimeException();

				//			List<Node> dcs = new LinkedList<> ();
				//			for (int d = 0; d < N; d ++)
				//				if (e_d.get(d) > 0) dcs.add(netPlan.getNode(d));
				//			if (dcs.size() != numOfDCs) throw new RuntimeException();

				//			if (!dcs.containsAll(previousDCPositionsThisCDN)) throw new RuntimeException();




				if ((bestSolution == null) || (bestSolution.getThird() > currentSolutionCost)) {
					bestSolution = Triple.of(r_uad, dcPositions, currentSolutionCost);
				}
				if (isFirstTime) flag = false;
			}
			flag = false;
		}while (flag);

		
		return Pair.of(bestSolution.getFirst(), bestSolution.getSecond());
	}
	
	
	
}