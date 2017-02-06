package com.net2plan.general;


import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;

import com.google.common.collect.Sets;
import com.net2plan.interfaces.networkDesign.Demand;
import com.net2plan.interfaces.networkDesign.IAlgorithm;
import com.net2plan.interfaces.networkDesign.Link;
import com.net2plan.interfaces.networkDesign.MulticastDemand;
import com.net2plan.interfaces.networkDesign.MulticastTree;
import com.net2plan.interfaces.networkDesign.NetPlan;
import com.net2plan.interfaces.networkDesign.Node;
import com.net2plan.libraries.GraphUtils;
import com.net2plan.libraries.TrafficMatrixGenerationModels;
import com.net2plan.utils.Constants.RoutingType;
import com.net2plan.utils.DoubleUtils;
import com.net2plan.utils.InputParameter;
import com.net2plan.utils.Pair;
import com.net2plan.utils.TimeTrace;
import com.net2plan.utils.Triple;

import cern.colt.matrix.tdouble.DoubleFactory1D;
import cern.colt.matrix.tdouble.DoubleFactory2D;
import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.jet.math.tdouble.DoubleFunctions;



/** This is a template to be used in the lab work, a starting point for the students to develop their programs
 * 
 */
public class TrafficTrendCDN implements IAlgorithm
{
	
	private InputParameter simYears = new InputParameter ("simYears", (int) 20 , "Simulation Years from 2016");
	private InputParameter G = new InputParameter ("G", (double) 1 , "Expansion Factor per CDN");
	private InputParameter sim = new InputParameter ("sim", (int) 1, "simulation index");
	private InputParameter isLocal = new InputParameter("isLocal", (boolean) false, "false if the execution is run in a server, true otherwise"); 
	private InputParameter A = new InputParameter("A", (int) 5, "Number of Apps");
	private InputParameter S = new InputParameter("S", (int) 5, "Number of Services");
	private InputParameter C = new InputParameter("C", (int) 5, "Number of CDNs");
	private InputParameter U = new InputParameter("U", (int) 100, "Number of Content Units");
	private InputParameter avNumReplicasPerContentUnit = new InputParameter("avNumReplicasPerContentUnit", (double) 2, "Average number of replicas of a content unit");
	private InputParameter rngSeed = new InputParameter("rngSeed", (long) 0, "The seed for the random numbers generator");
	private InputParameter debugMode = new InputParameter("debug Mode", (boolean) true, "Set up true if debug mode");

	private TimeTrace stat_averageRTTPerService_s = new TimeTrace();
	private TimeTrace stat_numberOfNewDCCreatedEachYear = new TimeTrace();
	private TimeTrace stat_sumTotalTrafficInLinksSummingOnlyTelcoTelco = new TimeTrace();
	private TimeTrace stat_sumTotalTrafficInLinksSummingOnlyDCToUserPerService_s = new TimeTrace();
	private TimeTrace stat_sumTotalTrafficInLinksSummingD2DPerService_s = new TimeTrace();
	private TimeTrace stat_sumTotalOfferedTrafficSummingOnlyDCToUserPerService_s = new TimeTrace();
	private TimeTrace stat_sumTotalOfferedTrafficSummingOnlyTelcoTelco = new TimeTrace();
	private TimeTrace stat_sumTotalOfferedTrafficSummingD2DPerService_s = new TimeTrace();

	//	private TimeTrace stat_sumTotalTrafficInTheLinks = new TimeTrace();
//	private TimeTrace stat_sumTotalOfferedTrafficOfDemandsEvenDemandsInSameNode = new TimeTrace();
	
	/** The method called by Net2Plan to run the algorithm (when the user presses the "Execute" button)
	 * @param netPlan The input network design. The developed algorithm should modify it: it is the way the new design is returned
	 * @param algorithmParameters Pair name-value for the current value of the input parameters
	 * @param net2planParameters Pair name-value for some general parameters of Net2Plan
	 * @return
	 */
	@SuppressWarnings("unchecked")
	@Override
	public String executeAlgorithm(NetPlan originalnetPlan, Map<String, String> algorithmParameters, Map<String, String> net2planParameters)
	{		
		long iniTime = System.nanoTime();
		originalnetPlan.removeAllDemands();
		
		/* Initialize all InputParameter objects defined in this object (this uses Java reflection) */
		InputParameter.initializeAllInputParameterFieldsOfObject(this, algorithmParameters);
		
		// Set the population of the nodes
//		NetPlan netPlan = originalnetPlan.copy();
		// Initial checks
		final int N = originalnetPlan.getNumberOfNodes();
		final int E = originalnetPlan.getNumberOfLinks();
		final int S = this.S.getInt();
		final int C = this.C.getInt();
		final int A = this.A.getInt();
		final int U = this.U.getInt();
		final int simYears = this.simYears.getInt();	
		final double H_D2C = 0.76*1000;
		final double H_TelcoTelco = 1000-H_D2C;
		final Random rand = new Random(this.rngSeed.getLong ());		

		double[] populationWeightVector = new double[N];
	
		final double CAGR_telcoTelco = 0.1;
		originalnetPlan.setRoutingType(RoutingType.SOURCE_ROUTING);	
		
		// Read services, applications and available CDNs
		TrafficTrendUtils appAndCDNInfo = new TrafficTrendUtils(C,S,A,U,rand);		
		appAndCDNInfo.generateServicePerApp();
		
		String path = null;
		if (!isLocal.getBoolean())
			path = "../libcplex1261.so";
		
		List<Triple<Double, Integer, int[]>> appInfo = appAndCDNInfo.getApps();
		List<List<Node>> nodesWithDCPerCDN_c = appAndCDNInfo.getInitialDCPlacementPerCDN(originalnetPlan , C , rand);
		List<Triple<Double,Double,Double>> servicesInfo = appAndCDNInfo.getServices();				
		
		int iniNumberDCs = 0;
		for(List<Node> dcsInThisCDN : nodesWithDCPerCDN_c)		
			iniNumberDCs += dcsInThisCDN.size();
			
		stat_numberOfNewDCCreatedEachYear.add(0,(double)iniNumberDCs);
		
		final DoubleMatrix1D populationPerNode = NetPlan.getAttributeValues(originalnetPlan.getNodes() , "population" , 0);
		final double totalPopulation = populationPerNode.zSum();
		for (int n = 0; n< N; n++)		
			populationWeightVector[n] = populationPerNode.get(n)/totalPopulation;			
		
		final DoubleMatrix1D linkCost = DoubleFactory1D.dense.make(E,1);
		final Pair<Double,Double> totalTrafficSummingTelcoTelcoYearZero = computeTelcoTelcoTotalTrafficInTheLinksYearZero (originalnetPlan , populationWeightVector , H_TelcoTelco , CAGR_telcoTelco);
		
		final Map<Pair<Node,Node>,List<List<Link>>> cpl = originalnetPlan.computeUnicastCandidatePathList (null , 1, -1, -1, -1, -1,-1, -1, null);
		final Map<List<Node>,Set<Link>> cplMulticast = new HashMap<> ();
		final DoubleMatrix2D rtt_n1n2 = DoubleFactory2D.dense.make(N,N);
		final DoubleMatrix2D numHops_n1n2 = DoubleFactory2D.dense.make(N,N);
		for (Pair<Node,Node> n1n2 : cpl.keySet())
		{
			final List<Link> seqLinks = cpl.get(n1n2).get(0);
			final double rtt = seqLinks.stream().mapToDouble(e->e.getPropagationDelayInMs()).sum();
			rtt_n1n2.set(n1n2.getFirst().getIndex() , n1n2.getSecond().getIndex() , rtt);
			numHops_n1n2.set(n1n2.getFirst().getIndex() , n1n2.getSecond().getIndex() , seqLinks.size());
		}
		List<List<List<List<Node>>>> replicaPlacements = appAndCDNInfo.computeReplicaPlacementsForAllCDNs (originalnetPlan , avNumReplicasPerContentUnit.getDouble() , rtt_n1n2 , populationWeightVector);
		
		for (int y = 0; y < simYears; y++)
		{
			
			double[] stat_sumTotalTrafficInLinksSummingOnlyDCToUserPerServiceThisYear_s = new double[S];
			double[] stat_sumTotalTrafficInLinksSummingD2DPerServiceThisYear_s = new double[S];
			double[] stat_sumTotalOfferedTrafficSummingOnlyDCToUserPerServiceThisYear_s = new double[S];
			double[] stat_sumTotalOfferedTrafficSummingD2DPerServiceThisYear_s = new double[S];
			double[] propagationTimeMultipliedByGbpsPerServiceThisYear_s = new double[S];

			stat_sumTotalOfferedTrafficSummingOnlyTelcoTelco.add(y,totalTrafficSummingTelcoTelcoYearZero.getFirst() * Math.pow(1 + CAGR_telcoTelco , y));
			stat_sumTotalTrafficInLinksSummingOnlyTelcoTelco.add(y,totalTrafficSummingTelcoTelcoYearZero.getSecond() * Math.pow(1 + CAGR_telcoTelco , y));
			
			// Compute traffic per application
			List<DoubleMatrix2D> traffMatrixThisYearIncludingSelfDemandsPerCDN_c = new ArrayList<DoubleMatrix2D>();
			for (int c = 0; c < C; c++)
				traffMatrixThisYearIncludingSelfDemandsPerCDN_c.add(c, DoubleFactory2D.sparse.make(N,N));			
			
			int appIndex = 0;
			for (Triple<Double,Integer, int[]> app : appInfo)
			{						
				final int appService = app.getSecond();																				// Type of service of the app
				final double cagr_service = servicesInfo.get(appService).getSecond();
				final double beta_s = servicesInfo.get(appService).getThird();
				final double appTraffic = H_D2C*app.getFirst()*Math.pow((1+cagr_service),y);
				final int[] indexesOfCDNsForThisApplication = app.getThird();												      	// List of CDNs where the app can carry the traffic				
				final double h_ac = appTraffic/(double) indexesOfCDNsForThisApplication.length;										// Traffic per each CDN for the app					
				
				for (int c=0; c < indexesOfCDNsForThisApplication.length; c ++)
				{
					for (int u=0; u < U; u++)					
					{
						// D2C Traffic Matrices
						List<Node> nodesWithDCWithAReplicaOfThisContentUnitInThisCDN = replicaPlacements.get(appIndex).get(c).get(u);
						final DoubleMatrix2D traffMatrixAppCDN = DoubleFactory2D.dense.make(N,N);
												
						for(Node n1 : originalnetPlan.getNodes())				
							for(Node n2 : nodesWithDCWithAReplicaOfThisContentUnitInThisCDN)						
								traffMatrixAppCDN.set(n1.getIndex(), n2.getIndex(), populationWeightVector[n1.getIndex()]*populationWeightVector[n2.getIndex()]);						
					
						final double sum_traffMatrixAppCDN = traffMatrixAppCDN.zSum();
						final DoubleMatrix2D normalizedTrafficMatrixAppCDN_nc = traffMatrixAppCDN.assign(DoubleFunctions.mult(h_ac*TrafficTrendUtils.zipfDitribution[u] / sum_traffMatrixAppCDN));
						traffMatrixThisYearIncludingSelfDemandsPerCDN_c.get(indexesOfCDNsForThisApplication[c]).assign(normalizedTrafficMatrixAppCDN_nc,DoubleFunctions.plus);
						/* sum the diagonal: this is offered traffic does not appear as demands */

						stat_sumTotalOfferedTrafficSummingOnlyDCToUserPerServiceThisYear_s[appService] += normalizedTrafficMatrixAppCDN_nc.zSum();

						/* The traffic goes 80% to the closest replica, and 20% spread among all (including closest replica) */
						for (Node userNode : originalnetPlan.getNodes ())
						{
							Node closestReplicaNode = null;
							int minNumHops = Integer.MAX_VALUE;
							for (Node n : nodesWithDCWithAReplicaOfThisContentUnitInThisCDN)
								if (numHops_n1n2.get(userNode.getIndex() , n.getIndex()) < minNumHops)
								{
									closestReplicaNode = n; 
									minNumHops = (int) numHops_n1n2.get(userNode.getIndex() , n.getIndex()); 
								}
							
							final double traffic = normalizedTrafficMatrixAppCDN_nc.get(userNode.getIndex() , closestReplicaNode.getIndex());

							/* 80% of traffic goes to closest replica */
							stat_sumTotalTrafficInLinksSummingOnlyDCToUserPerServiceThisYear_s[appService] += 0.8 * traffic * numHops_n1n2.get(userNode.getIndex() , closestReplicaNode.getIndex());
							propagationTimeMultipliedByGbpsPerServiceThisYear_s[appService] += 0.8 * traffic * rtt_n1n2.get(userNode.getIndex() , closestReplicaNode.getIndex());
							
							/* 20% of traffic is spread randomly among all the DCs in the CDN */
							final int numDCsWithReplicas = nodesWithDCWithAReplicaOfThisContentUnitInThisCDN.size()-1;

							for (Node n : nodesWithDCWithAReplicaOfThisContentUnitInThisCDN)
							{
								if (n.getIndex() != closestReplicaNode.getIndex())
								{
									stat_sumTotalTrafficInLinksSummingOnlyDCToUserPerServiceThisYear_s[appService] += (0.2 / numDCsWithReplicas)  * traffic * numHops_n1n2.get(userNode.getIndex() , n.getIndex());
									propagationTimeMultipliedByGbpsPerServiceThisYear_s[appService] += (0.2 / numDCsWithReplicas) * traffic * rtt_n1n2.get(userNode.getIndex() , n.getIndex());
								}
							}
						}

						// D2D Traffic Matrices (InterCDN Traffic)		
						final double Nmas = TrafficTrendUtils.zipfDitribution[u];
						final double replicaTraffic = beta_s*h_ac*Nmas;	
	
						final Node originNodeForMulticast = nodesWithDCWithAReplicaOfThisContentUnitInThisCDN.iterator().next();
						final Set<Node> destinationNodes = new HashSet<> (nodesWithDCWithAReplicaOfThisContentUnitInThisCDN);
						destinationNodes.remove(originNodeForMulticast);
						if(!destinationNodes.isEmpty())
						{
							Set<Link> multLinks = cplMulticast.get(nodesWithDCWithAReplicaOfThisContentUnitInThisCDN);
							if (multLinks == null)
							{
								multLinks = GraphUtils.getMinimumCostMulticastTree(originalnetPlan.getLinks(), originalnetPlan.getMatrixNodeLinkOutgoingIncidence(), originalnetPlan.getMatrixNodeLinkIncomingIncidence(),
									linkCost, originNodeForMulticast, destinationNodes, -1, -1, -1, -1, "cplex", path , 10);
								cplMulticast.put(nodesWithDCWithAReplicaOfThisContentUnitInThisCDN , multLinks);
							}
							stat_sumTotalOfferedTrafficSummingD2DPerServiceThisYear_s[appService] += replicaTraffic;
							stat_sumTotalTrafficInLinksSummingD2DPerServiceThisYear_s[appService] += replicaTraffic * multLinks.size();
						}
					}		
				}
				appIndex++;
			}
					
			for(int s=0; s<S; s++)			
				propagationTimeMultipliedByGbpsPerServiceThisYear_s[s] = propagationTimeMultipliedByGbpsPerServiceThisYear_s[s]/stat_sumTotalOfferedTrafficSummingOnlyDCToUserPerServiceThisYear_s[s];			
			
//			stat_sumTotalTrafficInTheLinks.add(y, netPlan.getVectorLinkCarriedTraffic().zSum());
//			stat_sumTotalOfferedTrafficOfDemandsEvenDemandsInSameNode.add(y, netPlan.getVectorDemandOfferedTraffic().zSum());
//			stat_sumTotalOfferedMulticastTraffic.add(y, netPlan.getVectorMulticastDemandOfferedTraffic().zSum());
			stat_averageRTTPerService_s.add(y, propagationTimeMultipliedByGbpsPerServiceThisYear_s);
			stat_sumTotalTrafficInLinksSummingOnlyDCToUserPerService_s.add(y, stat_sumTotalTrafficInLinksSummingOnlyDCToUserPerServiceThisYear_s);	
			stat_sumTotalOfferedTrafficSummingOnlyDCToUserPerService_s.add(y, stat_sumTotalOfferedTrafficSummingOnlyDCToUserPerServiceThisYear_s);
			stat_sumTotalTrafficInLinksSummingD2DPerService_s.add(y, stat_sumTotalTrafficInLinksSummingD2DPerServiceThisYear_s);	
			stat_sumTotalOfferedTrafficSummingD2DPerService_s.add(y, stat_sumTotalOfferedTrafficSummingD2DPerServiceThisYear_s);

			int newDC = 0;	
			// Update Data Center Locations
			if(G.getDouble() > 0)	
			for(int c = 0; c < C; c++)			
				newDC += appAndCDNInfo.addDcIntoCDN(originalnetPlan, c, G.getDouble(),traffMatrixThisYearIncludingSelfDemandsPerCDN_c.get(c));
			
			stat_numberOfNewDCCreatedEachYear.add(y+1,newDC);
			replicaPlacements = appAndCDNInfo.computeReplicaPlacementsForAllCDNs (originalnetPlan , avNumReplicasPerContentUnit.getDouble() , rtt_n1n2 , populationWeightVector);
		
			if(debugMode.getBoolean())
			{
				double timenow = (System.nanoTime()-iniTime)*1e-9; 
				System.out.println("---------- Year " + Integer.toString(y) + " ----------");
				System.out.println(" Year "+y+ " in " + timenow + " seconds.");
				System.out.println("* Non CDN Traffic: ");
				System.out.println("  - Offered: " + totalTrafficSummingTelcoTelcoYearZero.getFirst() * Math.pow(1 + CAGR_telcoTelco , y));
				System.out.println("  - In Links: " + totalTrafficSummingTelcoTelcoYearZero.getSecond() * Math.pow(1 + CAGR_telcoTelco , y));
				System.out.println("* DC to User Traffic: ");
				System.out.println("  - Offered: " + DoubleUtils.sum(stat_sumTotalOfferedTrafficSummingOnlyDCToUserPerServiceThisYear_s));
				System.out.println("  - In links: " +  DoubleUtils.sum(stat_sumTotalTrafficInLinksSummingOnlyDCToUserPerServiceThisYear_s));
				System.out.println("* DC to DC traffic: " );
				System.out.println("  - Offered: " +  DoubleUtils.sum(stat_sumTotalOfferedTrafficSummingD2DPerServiceThisYear_s));
				System.out.println("  - In links: " +  DoubleUtils.sum(stat_sumTotalTrafficInLinksSummingD2DPerServiceThisYear_s));
				System.out.println("* RTT : " + DoubleUtils.sum(propagationTimeMultipliedByGbpsPerServiceThisYear_s));
				System.out.println("* Num New DCs: " + newDC);
			}
			
		}
					
		String root;
		if (isLocal.getBoolean()) root = "C:/Users/Javi/OneDrive/Projects/Proyecto Trend CDN/Results/";
		else root = "../trendTraffic/Results/";
		
		String gString = Double.toString(G.getDouble());
		String simString = Integer.toString(sim.getInt());
		
		stat_averageRTTPerService_s.printToFile(new File(root+"avRtt_"+gString+"_"+simString+".txt"));	
		stat_sumTotalTrafficInLinksSummingOnlyTelcoTelco.printToFile(new File(root+"carriedNotCDNTraffic_"+gString+"_"+simString+".txt"));		
		stat_sumTotalOfferedTrafficSummingOnlyTelcoTelco.printToFile(new File(root+"offeredNotCDNTraffic_"+gString+"_"+simString+".txt"));
		stat_sumTotalTrafficInLinksSummingOnlyDCToUserPerService_s.printToFile(new File(root+"carriedD2CTrafficPerService_"+gString+"_"+simString+".txt"));
		stat_sumTotalOfferedTrafficSummingOnlyDCToUserPerService_s.printToFile(new File(root+"offeredD2CTrafficPerService_"+gString+"_"+simString+".txt"));
		stat_sumTotalOfferedTrafficSummingD2DPerService_s.printToFile(new File(root+"multicast_offered"+gString+"_"+simString+".txt"));	
		stat_sumTotalTrafficInLinksSummingD2DPerService_s.printToFile(new File(root+"multicast_inlinks"+gString+"_"+simString+".txt"));
		stat_numberOfNewDCCreatedEachYear.printToFile(new File(root+"numberOfNewDC"+gString+"_"+simString+".txt"));
		
		return "Ok!"; // this is the message that will be shown in the screen at the end of the algorithm
	}

	/** Returns a description message that will be shown in the graphical user interface
	 */
	@Override
	public String getDescription()
	{
		return "Here you should return the algorithm description to be printed by Net2Plan graphical user interface";
	}

	
	/** Returns the list of input parameters of the algorithm. For each parameter, you should return a Triple with its name, default value and a description
	 * @return
	 */
	@Override
	public List<Triple<String, String, String>> getParameters()	
	{
		return InputParameter.getInformationAllInputParameterFieldsOfObject(this);
	}

	private Pair<Double,Double> computeTelcoTelcoTotalTrafficInTheLinksYearZero (final NetPlan np , final double [] popuWeightVector , final double H_TelcoTelco , final double CAGR_telcoTelco)
	{
		final NetPlan netPlan = np.copy();
		final int N = netPlan.getNumberOfNodes();
		netPlan.removeAllDemands();

		DoubleMatrix2D traffMatrixTelcoTelco = DoubleFactory2D.dense.make(N,N);	
		for(int n1 = 0; n1 < N; n1++)				
			for(int n2 = 0; n2 < N; n2++)						
				if(n1 != n2)				
					traffMatrixTelcoTelco.set(n1, n2, popuWeightVector[n1]*popuWeightVector[n2]);					

		final double telecoTelcoTraffic = H_TelcoTelco;
		final DoubleMatrix2D initialTraffMatrixTelcoTelco = TrafficMatrixGenerationModels.normalizationPattern_totalTraffic(traffMatrixTelcoTelco, telecoTelcoTraffic);
		final List<Demand> notCDNDmands = netPlan.addDemandsFromTrafficMatrix(initialTraffMatrixTelcoTelco);			
		
		for (Demand d: notCDNDmands)
		{
			final List<Link> seqLinks = GraphUtils.getShortestPath(netPlan.getNodes() ,netPlan.getLinks(), d.getIngressNode(), d.getEgressNode(), null);
			netPlan.addRoute(d, d.getOfferedTraffic(), d.getOfferedTraffic(), seqLinks, null);
		}
		return Pair.of(netPlan.getVectorDemandOfferedTraffic().zSum() , netPlan.getVectorLinkCarriedTraffic().zSum());
	}
	
}
