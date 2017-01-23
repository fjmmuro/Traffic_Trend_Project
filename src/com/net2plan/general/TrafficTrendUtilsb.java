package com.net2plan.general;


import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Random;

import org.apache.commons.collections15.map.HashedMap;

import com.jom.DoubleMatrixND;
import com.jom.OptimizationProblem;
import com.net2plan.interfaces.networkDesign.Demand;
import com.net2plan.interfaces.networkDesign.Link;
import com.net2plan.interfaces.networkDesign.Net2PlanException;
import com.net2plan.interfaces.networkDesign.NetPlan;
import com.net2plan.interfaces.networkDesign.Node;
import com.net2plan.utils.Pair;
import com.net2plan.utils.Quadruple;
import com.net2plan.utils.RandomUtils;
import com.net2plan.utils.Triple;

import cern.colt.matrix.tdouble.DoubleFactory2D;
import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;


public class TrafficTrendUtilsb 
{
	
	private int C = 5;
	private int S = 5;
	private int A = 5;
	private int U = 100;
	
	
	private double[] x_s = {0.47,0.25,0.19,0.08,0.01}; 		// Hay que cambiarlo por los valores "reales del cisco vni" no se si normalizarlo a 1 o ponerlo como proporcion del total de trafico
	private double[] cagr = {0.31,0.31,0.18,0,0.47};  		// Cagr for each service regarging cisco vni
	private double[] beta = {0.1,0.5,1,0.9,0.1};			// Deberian valores menores
	private double[] N_u = new double[U];
	
	double[] lastInitialTraffic = new double[C];
	private Random rand = new Random();
	
	DoubleMatrix2D xas = DoubleFactory2D.dense.make(A,S);
	List<List<Integer>> cdnIndexes = new ArrayList<List<Integer>>();
	List<List<List<Integer>>> replicaPlacements = new ArrayList<List<List<Integer>>>();
	private List<Triple<Double, Double,Double>> services = new ArrayList<Triple<Double,Double,Double>>();
	
	public TrafficTrendUtilsb ()
	{		
		for(int c = 0; c < C; c++)
			lastInitialTraffic[c] = 0;
	}
	
	// This method set the node population
		public void setNodePopulation(NetPlan netPlan)
		{
			for(Node n : netPlan.getNodes())
			{
				String state = n.getAttribute("state");
				String population = "population";
				switch (state) 
				{
					case "Texas":
						n.setAttribute(population, "4542417");
						break;
					case "New York":
						n.setAttribute(population, "3299298");
						break;
					case "New Mexico":
						n.setAttribute(population, "2085109");
						break;
					case "Georgia":
						n.setAttribute(population, "10214860");
						break;
					case "Maryland":
						n.setAttribute(population, "6006401");
						break;
					case "Louisiana":
						n.setAttribute(population, "2335362");
						break;
					case "Montana":
						n.setAttribute(population, "1032949");
						break;
					case "Alabama":
						n.setAttribute(population, "4858979");
						break;
					case "Massachusetts":
						n.setAttribute(population, "6794422");
						break;
					case "North Dakota":
						n.setAttribute(population, "756927");
						break;
					case "South Carolina":
						n.setAttribute(population, "4896146");
						break;
					case "North Carolina":
						n.setAttribute(population, "3347600");
						break;
					case "Illinois":
						n.setAttribute(population, "6429997");
						break;
					case "Ohio":
						n.setAttribute(population, "2903355");
						break;
					case "Michigan":
						n.setAttribute(population, "9922576");
						break;
					case "Connecticut":
						n.setAttribute(population, "3590886");
						break;
					case "Florida":
						n.setAttribute(population, "3378545");
						break;
					case "Missouri":
						n.setAttribute(population, "3041836");
						break;
					case "Nevada":
						n.setAttribute(population, "2890845");
						break;
					case "Arkansas":
						n.setAttribute(population, "2978204");
						break;
					case "Kentucky":
						n.setAttribute(population, "4425092");
						break;
					case "Tennessee":
						n.setAttribute(population, "6600299");
						break;
					case "Wisconsin":
						n.setAttribute(population, "5771337");
						break;
					case "Minnesota":
						n.setAttribute(population, "5489594");
						break;
					case "New Jersey":
						n.setAttribute(population, "8958013");
						break;
					case "Virginia":
						n.setAttribute(population, "2794331");
						break;
					case "Oklahoma":
						n.setAttribute(population, "1955669");
						break;
					case "Nebraska":
						n.setAttribute(population, "1896190");
						break;
					case "Pennsylvania":
						n.setAttribute(population, "4267501");
						break;
					case "Arizona":
						n.setAttribute(population, "3414032");
						break;
					case "Oregon":
						n.setAttribute(population, "4028977");
						break;
					case "Rhode Island":
						n.setAttribute(population, "1056298");
						break;
					case "Utah":
						n.setAttribute(population, "2995919");
						break;
					case "Washington":
						n.setAttribute(population, "3585175");
						break;
					case "Delaware":
						n.setAttribute(population, "945934");
						break;
					case "California":
						n.setAttribute(population, "4893102");
						break;
					case "Colorado":
						n.setAttribute(population, "5456574");
						break;
				}
			}	
		}
	
	// Function to generate the CDN randomly (list<indexes>)
	public List<List<Integer>> getCDNs(int N)
	{		
		for(int c = 0; c < C; c++)
		{
			//int numNodes = RandomUtils.random(2, 4);	
			int numNodes = 3;	
			List<Integer> indexesThisCDN = new ArrayList<Integer>();
			for (int n = 0; n < numNodes; n++)
			{				 
				int candidateIndex = RandomUtils.random(0, N-1);
				while(indexesThisCDN.contains(candidateIndex))
					candidateIndex = RandomUtils.random(0, N-1);				
				indexesThisCDN.add(candidateIndex);
			}
			this.cdnIndexes.add(indexesThisCDN);
		}		
		return this.cdnIndexes;
	}
	
	// Function to generate the services randomly (x_s,cagr)
	public List<Triple<Double, Double,Double>> getServices()
	{
		for(int s = 0; s< S; s++)		
			this.services.add(new Triple<Double, Double,Double>(this.x_s[s], this.cagr[s], this.beta[s], false));
				
		return this.services;
	}
	
	
	// Function to generate the applications randomly (x_as,service,cdns)
	public List<Quadruple<Double,Integer,Integer,double[]>> getApps(int U, String solver, String solverName)
	{
		List<Quadruple<Double, Integer, Integer, double[]>> apps = new ArrayList<Quadruple<Double,Integer,Integer,double[]>>();
		List<Integer> shuffleCDNs = new ArrayList<Integer>();
		for(int c = 0; c < C; c++)
			shuffleCDNs.add(c);
		
		Pair<List<Integer>,int[]> services = generateServicePerApp();
		List<Integer> indexAppInService = services.getFirst();
		int[] numAppPerService = services.getSecond();
//		double[] x_a = generateAppTrafficPerService(numAppPerService,indexAppInService);		
		
		for(int a = 0; a < A; a++)
		{
			double[] N_u = computeNumberOfAccesses(U);
			int numCDNs = RandomUtils.random(1, 3);
			Collections.shuffle(shuffleCDNs);
			int[] appCDNs = new int[numCDNs];
			
			for (int c = 0; c < numCDNs; c++)
				appCDNs[c] = shuffleCDNs.get(c);

			apps.add(Quadruple.of(x_s[a], a, a,N_u));
		}		
		
		return apps;
	}	
	
	public double[] generateAppTrafficPerService(int[] numApp, List<Integer> indexAppInService)
	{
		double x_as[] = new double[A];		
		double[] totalX_a = new double[S];
		int[] numAppAux = new int[S];
		
		for(int a = 0; a < A; a++)
		{
			int s = indexAppInService.get(a);
			if (numAppAux[s] == numApp[s]-1)			
				x_as[a] = 1- totalX_a[s];			
			else
			{
				double x_a = rand.nextDouble();
				while((x_a + totalX_a[s] > 1))
					x_a = rand.nextDouble();				
				numAppAux[s]++;
				totalX_a[s] += x_a;
				x_as[a] = x_a;
			}			
		}
		
		return x_as;
	}
	
	public Pair<List<Integer>,int[]> generateServicePerApp()
	{
		int[] numAppPerService = new int[S];
		List<Integer> indexAppInService = new ArrayList<Integer>();
		boolean isAllServicesGenerated = false;		
		
		while(!isAllServicesGenerated)
			for(int a = 0; a < A; a++)
			{			
				int service = RandomUtils.random(0, S-1);
				indexAppInService.add(service);
				numAppPerService[service] ++;
	
				for (int s = 0; s < S; s ++)	
				{
					if(!indexAppInService.contains(s))
						break;
					else if (indexAppInService.contains(S-1))
						isAllServicesGenerated = true;
					else
						break;							
				}		
			}
		
		return Pair.of(indexAppInService, numAppPerService);		
		
	}
	
	public int addDcIntoCDN(NetPlan netPlan, int c, double G, DoubleMatrix2D trafficMatrix)
	{
		List<Integer> currentIndexes = this.cdnIndexes.get(c);

		int numberOfNewDCs = 0;
		
		if (this.lastInitialTraffic[c] == 0) this.lastInitialTraffic[c] =  trafficMatrix.zSum();
		double intialCDNsTraffic = this.lastInitialTraffic[c];
		double currentCDNTraffic = trafficMatrix.zSum();
		int N_dc = (int) (G*(currentCDNTraffic-intialCDNsTraffic)/intialCDNsTraffic);	
			
		int N = netPlan.getNumberOfNodes();
		
		if(N_dc >= 1 && currentIndexes.size()+N_dc <= N)
		{
			for (int n1 = 0; n1 < N_dc; n1++)
			{
				double bestTraffic = 0;
				int bestIndex = -1;
				double nodeTraffic = 0;
				
				if (rand.nextDouble() < 0.8)					
				{
					for (int n = 0; n < N; n++)
					{
						if (!currentIndexes.contains(n))
						{								
							Node node = netPlan.getNode(n);							
							nodeTraffic = node.getIngressCarriedTraffic()+node.getEgressCarriedTraffic();
							
							if (nodeTraffic > bestTraffic)
							{
								bestTraffic = nodeTraffic;
								bestIndex = n;
							}	
						}
					}
				}
				else
				{
					int newCandidate = RandomUtils.random(0, N-1);
					while (currentIndexes.contains(newCandidate))
						newCandidate = RandomUtils.random(0, N-1);					
					bestIndex = newCandidate;					
				}	
				if(bestIndex != -1)
				{
					currentIndexes.add(bestIndex);	
					updateReplicaPlacements(c,bestIndex);
					numberOfNewDCs++;
				}
			}
			this.lastInitialTraffic[c] = currentCDNTraffic;
		}
		return numberOfNewDCs;
		
	}
	
	public double[] computeNumberOfAccesses(int U)
	{		
		double[] x_u = new double[U];
		double total = 0;
		
		for(int i=0; i < U; i++)
		{
			x_u[i] = (double) U/(i+1);
			total += x_u[i];
		}
		for (int i = 0; i < U; i++)
			N_u[i] = x_u[i]/total;	
		
		return this.N_u;
	}
	
	public List<List<List<Integer>>> computeInitialReplicaPlacements (List<List<Integer>> cdnIndexesNodes)
	{
		for(List<Integer> cdn : cdnIndexesNodes)
		{
			List<List<Integer>> replicaPlacesThisCDN = new ArrayList<List<Integer>>();
			int cdnSize = cdn.size();
			
			for(int i = 0; i < U; i++)			
			{
				int numReplicas = RandomUtils.random(2, cdnSize);
				List<Integer> cdnCopy = cdn;
				Collections.shuffle(cdnCopy);
				List<Integer> indexesContentUnit = new ArrayList<Integer>();
				for(int index = 0; index < numReplicas; index++)
					indexesContentUnit.add(cdnCopy.get(index));
				replicaPlacesThisCDN.add(indexesContentUnit);
			}		
			this.replicaPlacements.add(replicaPlacesThisCDN);			
		}
			
		return this.replicaPlacements;
	}
	
	public void updateReplicaPlacements(int c, int bestIndex)
	{
		for (int i = 0; i < U; i++)
			if (rand.nextDouble() < 200*this.N_u[i])
				this.replicaPlacements.get(c).get(i).add(bestIndex);			
	}
	
	public List<List<List<Integer>>> checkReplicaPlacements()
	{
		return this.replicaPlacements;
	}
	
}
		