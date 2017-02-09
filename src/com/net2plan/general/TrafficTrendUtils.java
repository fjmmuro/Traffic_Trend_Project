package com.net2plan.general;

/* 
 * Hypothesis:
// * 1) average number of replicas of CU is proportional to its popularity. avNumberReplicasInCDN (u) = beta * zipf(u)
// * avNumberReplicasInCDN = sum (u) beta * zipf(u) = beta
 * 2) probability of a DC of having a replica of CU is proportional to its popularity: probHavingReplicaofCUInsideDC = alpha * zipf(u) 
 * 
 * avNumberReplicasInCDN (u) = numDCssOfCDN * probHavingReplicaofCUInsideDC (u) =>
 * => sum(u) avNumberReplicasInCDN = numDCsOfCDN sum (u) probHavingReplicaofCUInsideDC (u)
 * => beta * sum(u) zipf(u) = numDCsOfCDN alpha sum (u) zipf (u)
 * => beta =  numDCsOfCDN * alpha
 * 
 * Normalize so that AvNumberOfReplicasOfContentUnitInADC = XXX
 * 
 * avNumberReplicasInCDN = 100 * XXX => beta = 100 * XXX => alpha = 100 * XXX /  numDCsOfCDN
 * 
 * */



import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Random;
import java.util.Set;

import com.google.common.collect.Sets;
import com.net2plan.interfaces.networkDesign.NetPlan;
import com.net2plan.interfaces.networkDesign.Node;
import com.net2plan.utils.RandomUtils;
import com.net2plan.utils.Triple;

import cern.colt.matrix.tdouble.DoubleMatrix2D;

public class TrafficTrendUtils 
{
	
	private int C,S,A,U;
	
	private double[] x_s = {0.47,0.25,0.19,0.08,0.01}; 		// Total traffic proportion for each service
	private double[] cagr = {0.31,0.31,0.18,0,0.47};  		// Cagr for each service regarding cisco vni
	private double[] beta = {0.1,0.5,1,0.9,0.1};			// Beta values for each service
	
	public static double [] zipfDitribution = null;
	
	double[] trafficInPreivousYearWhenADCWasCreated;
	Boolean[] isModifiedCDN;
	private Random rand;
	
	private List<List<Node>> cdnNodes_c = new ArrayList<>();
	private List<Triple<Double, Double,Double>> services = new ArrayList<>();
	private List<Triple<Double, Integer, int[]>> apps = new ArrayList<>();
	private List<List<List<List<Node>>>> lastYearReplicaPlacements = new ArrayList<>();	
	
	public TrafficTrendUtils (int C,int S,int A,int U,Random rand)
	{		
		this.rand = rand;
		this.C = C;
		this.S = S;
		this.A = A;
		this.U = U;
		
		this.trafficInPreivousYearWhenADCWasCreated = new double[C];
		this.isModifiedCDN = new Boolean[C];
		
		for(int c = 0; c < C; c++)
		{
			trafficInPreivousYearWhenADCWasCreated[c] = 0;
			isModifiedCDN[c] = false;
		}
		
		TrafficTrendUtils.updateZipfDistributionValues(U);
	}
	
	// Function to generate the CDN randomly (list<indexes>)
	public List<List<Node>> getInitialDCPlacementPerCDN(NetPlan np , final int C , final Random rand)
	{		
		final List<List<Node>> res  = new ArrayList<> ();
		for(int c = 0; c < C; c++)
		{
			final int numNodes = 3;	
			final List<Node> shuffledNodes = new ArrayList<> (np.getNodes());
			Collections.shuffle(shuffledNodes , rand);
			res.add(new ArrayList<> (shuffledNodes.subList(0 , numNodes)));
			this.isModifiedCDN[c] = true;
		}		
		this.cdnNodes_c = res;
		
		return this.cdnNodes_c;
	}
	
	// Function to generate the services randomly (x_s,cagr)
	public List<Triple<Double, Double,Double>> getServices()
	{
		for(int s = 0; s< S; s++)		
			this.services.add(new Triple<Double, Double,Double>(this.x_s[s], this.cagr[s], this.beta[s], false));
				
		return this.services;
	}
	
	
	// Function to generate the applications randomly (x_as,service,cdns)
	public List<Triple<Double,Integer,int[]>> getApps()
	{
		List<Integer> appServiceIndex = generateServicePerApp();
		List<Integer> shuffleCDNs = new ArrayList<Integer>();
		for(int c = 0; c < C; c++)
			shuffleCDNs.add(c);
		
		for(int a = 0; a < A; a++)		
		{			
			int appService = appServiceIndex.get(a);
			int numberofAppsPerThisService = Collections.frequency(appServiceIndex, appService);
			int numCDNsForThisApplication = RandomUtils.random(1, 3);
			Collections.shuffle(shuffleCDNs);
			int[] indexesOfCDNsOfThisApplication = new int[numCDNsForThisApplication];
			double x_a = x_s[appService]/(double) numberofAppsPerThisService;
			
			for (int c = 0; c < numCDNsForThisApplication; c++)
				indexesOfCDNsOfThisApplication[c] = shuffleCDNs.get(c);

			this.apps.add(Triple.of(x_a, appService, indexesOfCDNsOfThisApplication));
		}		
		
		return this.apps;
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
	
	public List<Integer> generateServicePerApp()
	{
		List<Integer> indexAppInService = new ArrayList<Integer>();	

		for(int a = 0; a < A; a++)			
			if (a < S)
				indexAppInService.add(a);
			else
				indexAppInService.add(RandomUtils.random(0, S-1));		
		
		return indexAppInService;				
	}
	
	public int addDcIntoCDN(NetPlan netPlan, int c, double G, DoubleMatrix2D trafficMatrix)
	{
		final List<Node> originalDCsInCDN = new ArrayList<Node> (this.cdnNodes_c.get(c)); 
		List<Node> currentDCsInCDN = this.cdnNodes_c.get(c);

		/* First time passes */
		if (this.trafficInPreivousYearWhenADCWasCreated[c] == 0) this.trafficInPreivousYearWhenADCWasCreated[c] =  trafficMatrix.zSum();

		double intialCDNsTraffic = this.trafficInPreivousYearWhenADCWasCreated[c];
		double currentCDNTraffic = trafficMatrix.zSum();
		
		final int numberOfNewDCsToCreate = (int) (G*(currentCDNTraffic-intialCDNsTraffic)/intialCDNsTraffic);			
		final int N = netPlan.getNumberOfNodes();
				
		if(numberOfNewDCsToCreate >= 1 && currentDCsInCDN.size() + numberOfNewDCsToCreate <= N)
		{
			for (int n1 = 0; n1 < numberOfNewDCsToCreate; n1++)
			{
				/* Selection of the placement for the new DC */
				final Set<Node> placementCandidates = Sets.difference(new HashSet<>(netPlan.getNodes()) , new HashSet<>(currentDCsInCDN));
				Node chosenNode = null;
				if (rand.nextDouble() < 0.8)					
				{
					/* Take the best option */
					double bestTraffic = -Double.MAX_VALUE;
					for (Node candidate : placementCandidates)
					{
						final double nodeTraffic = trafficMatrix.viewColumn(candidate.getIndex()).zSum() + trafficMatrix.viewRow(candidate.getIndex()).zSum();   
						if (nodeTraffic > bestTraffic)
						{
							bestTraffic = nodeTraffic;
							chosenNode = candidate;
						}	
					}
				}
				else
					chosenNode = new ArrayList<Node> (placementCandidates).get(rand.nextInt(placementCandidates.size()));
				
				if(chosenNode == null) throw new RuntimeException("Chosen Node = null");
				this.cdnNodes_c.get(c).add(chosenNode);
			}
			this.trafficInPreivousYearWhenADCWasCreated[c] = currentCDNTraffic;
			this.isModifiedCDN[c] = true;
		}		
		final int numberOfNewDCs = this.cdnNodes_c.get(c).size() - originalDCsInCDN.size() ;
		
		
		return  numberOfNewDCs;
	}
	
	/** For making the Zipf distribution 
	 * @param U
	 * @return
	 */
	public static void updateZipfDistributionValues (int U)
	{		
		zipfDitribution = new double [U];
		double[] x_u = new double[U];
		double total = 0;
		
		for(int i=0; i < U; i++)
		{
			x_u[i] = (double) U/(i+1);
			total += x_u[i];
		}
		for (int i = 0; i < U; i++)
			zipfDitribution [i] = x_u[i]/total;	
	}
	
	public List<List<List<List<Node>>>> computeReplicaPlacementsForAllCDNs (NetPlan np , double averageNumberOfReplicasPerCU , DoubleMatrix2D rtt_n1n2 , double [] population_n)
	{
		
		List<List<List<List<Node>>>> replicaPlacements_acu = new ArrayList<>();
		int a = 0;
		for(Triple<Double, Integer, int[]> apps : this.apps)
		{
			List<List<List<Node>>> replicaPlacementsThisApp = new ArrayList<>();
			int[] cdnsThisApp = apps.getThird();
			for(int c = 0; c < cdnsThisApp.length; c++)
			{			
				
				if( (lastYearReplicaPlacements.isEmpty()) || (this.isModifiedCDN[c]) )
				{
					List<Node> cdnDCsThisYear = this.cdnNodes_c.get(cdnsThisApp[c]);
					List<List<Node>> replicaPlacesThisCDNAllCUs = new ArrayList<>();
					final int cdnNumDCs = cdnDCsThisYear.size();
					final int maximumNumberReplicasInEachDCEachApp = (int) Math.ceil(averageNumberOfReplicasPerCU * U) ;
					
//					double iniTime = (double) System.nanoTime()*1e-9;
					DoubleMatrix2D replicasPlacementsInThisCDN = ReplicaPlacement.placeReplicas(np, cdnDCsThisYear, rtt_n1n2, maximumNumberReplicasInEachDCEachApp, zipfDitribution,population_n , U); 
//					double endTIme = (double) System.nanoTime()*1e-9;
//					double ilpTime = endTIme-iniTime;
//					System.out.println("ILP Run time: " + ilpTime);
					
					
					for(int u = 0; u < U; u++)			
					{
						List<Node> dcsThisCU = new ArrayList<>();
						for(int d = 0; d < cdnNumDCs; d++ )
							if(replicasPlacementsInThisCDN.get(u, d) == 1)													
								dcsThisCU.add(cdnDCsThisYear.get(d));	
						replicaPlacesThisCDNAllCUs.add(dcsThisCU);
					}		
					replicaPlacementsThisApp.add(replicaPlacesThisCDNAllCUs);		
				}
				else				
					replicaPlacementsThisApp.add(lastYearReplicaPlacements.get(a).get(c));								
			}
			replicaPlacements_acu.add(replicaPlacementsThisApp);
			a++;
		}		
		for(int c = 0; c < C; c++)
			this.isModifiedCDN[c] = false;
		this.lastYearReplicaPlacements = replicaPlacements_acu;
		
		return replicaPlacements_acu;
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
		
		
		
	
}
		