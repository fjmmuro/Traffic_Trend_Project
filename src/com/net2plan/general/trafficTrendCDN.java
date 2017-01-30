package com.net2plan.general;

import java.io.File;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;

import com.net2plan.interfaces.networkDesign.Demand;
import com.net2plan.interfaces.networkDesign.IAlgorithm;
import com.net2plan.interfaces.networkDesign.Link;
import com.net2plan.interfaces.networkDesign.MulticastDemand;
import com.net2plan.interfaces.networkDesign.MulticastTree;
import com.net2plan.interfaces.networkDesign.NetPlan;
import com.net2plan.interfaces.networkDesign.Node;
import com.net2plan.interfaces.networkDesign.Route;
import com.net2plan.utils.InputParameter;
import com.net2plan.utils.TimeTrace;
import com.net2plan.utils.Triple;
import com.net2plan.utils.Constants.RoutingType;
import com.net2plan.libraries.GraphUtils;
import com.net2plan.libraries.TrafficMatrixGenerationModels;

import cern.colt.matrix.tdouble.DoubleFactory1D;
import cern.colt.matrix.tdouble.DoubleFactory2D;
import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.jet.math.tdouble.DoubleFunctions;


/** This is a template to be used in the lab work, a starting point for the students to develop their programs
 * 
 */
public class trafficTrendCDN implements IAlgorithm
{
	
	private InputParameter simYears = new InputParameter ("simYears", (int) 20 , "Simulation Years from 2016");
	private InputParameter G = new InputParameter ("G", (double) 1 , "Expansion Factor per CDN");
	private InputParameter sim = new InputParameter ("sim", (int) 1, "simulation index");
	private InputParameter isLocal = new InputParameter("isLocal", (boolean) false, "false if the execution is run in a server, true otherwise"); 
	private InputParameter A = new InputParameter("A", (int) 5, "Number of Apps");
	private InputParameter S = new InputParameter("S", (int) 5, "Number of Services");
	private InputParameter C = new InputParameter("C", (int) 5, "Number of CDNs");
	private InputParameter U = new InputParameter("U", (int) 100, "Number of Content Units");

	private TimeTrace stat_avRTT = new TimeTrace();
	private TimeTrace stat_carriedNotCDNTraffic = new TimeTrace();
	private TimeTrace stat_carriedD2CTrafficPerService = new TimeTrace();
	private TimeTrace stat_offeredD2CTrafficPerService = new TimeTrace();
	private TimeTrace stat_numberOfNewDC = new TimeTrace();
	private TimeTrace stat_multiCastTraffic = new TimeTrace();
	private TimeTrace stat_carried = new TimeTrace();
	private TimeTrace stat_offered = new TimeTrace();
	
	/** The method called by Net2Plan to run the algorithm (when the user presses the "Execute" button)
	 * @param netPlan The input network design. The developed algorithm should modify it: it is the way the new design is returned
	 * @param algorithmParameters Pair name-value for the current value of the input parameters
	 * @param net2planParameters Pair name-value for some general parameters of Net2Plan
	 * @return
	 */
	@Override
	public String executeAlgorithm(NetPlan originalnetPlan, Map<String, String> algorithmParameters, Map<String, String> net2planParameters)
	{		
//		long iniTime = System.nanoTime();
		originalnetPlan.removeAllDemands();
		
		/* Initialize all InputParameter objects defined in this object (this uses Java reflection) */
		InputParameter.initializeAllInputParameterFieldsOfObject(this, algorithmParameters);
		
		// Set the population of the nodes
		NetPlan netPlan = originalnetPlan.copy();
		// Initial checks
		final int N = netPlan.getNumberOfNodes();
		final int S = this.S.getInt();
		final int C = this.C.getInt();
		final int A = this.A.getInt();
		final int U = this.U.getInt();
		final int simYears = this.simYears.getInt();	
		final double H_D2C = 0.76*1000;
		final double H_TelcoTelco = 1000-H_D2C;
		final Random rand = new Random();		

		int[] populationVector = new int[N];
		double[] popuWeightVector = new double[N];
		int totalPopulation = 0;
		int[] levelVector = new int[N];						
	
		double CAGR_telcoTelco = 0.1;
		double telecoTelcoTraffic = H_TelcoTelco;
		double Nmas, replicaTraffic;
		double beta_s;
		double appTraffic;    					// Total traffic of the application app
		double h_ac;	
		double cagr_service;

		int originNodeReplicaIndex;		
		int numberOfCDNs;													// Number of different CDN available in this app
		int appService;														// Type of service of the app		
		
		List<Link> seqLinks = new ArrayList<Link>();
		Set<Link> multLinks = new HashSet<Link>();
		
		MulticastDemand multicastDemand;
		MulticastTree multicastTree;
		
		netPlan.setRoutingType(RoutingType.SOURCE_ROUTING);	
		
		// Read services, applications and available CDNs
		TrafficTrendUtils appAndCDNInfo = new TrafficTrendUtils(C,S,A,U, rand);		
//		appAndCDNInfo.setNodePopulation(netPlan);	
		appAndCDNInfo.generateServicePerApp();
		
		String path = null;
		if (!isLocal.getBoolean())
			path = "../libcplex1261.so";
		
		List<Triple<Double, Integer, int[]>> appInfo = appAndCDNInfo.getApps();
		List<List<Integer>> cdnIndexesNodes = appAndCDNInfo.getCDNs(N);
		List<Triple<Double,Double,Double>> servicesInfo = appAndCDNInfo.getServices();				
		
		int iniNumberDCs = 0;
		for(List<Integer> cdn : cdnIndexesNodes)		
			iniNumberDCs += cdn.size();
			
		stat_numberOfNewDC.add(0,(double)iniNumberDCs/C);
		
		for (Node n: netPlan.getNodes())
		{
			levelVector[n.getIndex()] = 1;
			populationVector[n.getIndex()] = Integer.parseInt(n.getAttribute("population"));
			totalPopulation += populationVector[n.getIndex()];
		}
				
		for (int n = 0; n< N; n++)		
			popuWeightVector[n] = (double)populationVector[n]/(double)totalPopulation;			
		
		final DoubleMatrix1D linkCost = DoubleFactory1D.dense.make(netPlan.getNumberOfLinks(),1);
		DoubleMatrix2D traffMatrixTelcoTelco = DoubleFactory2D.dense.make(N,N);	
		for(int n1 = 0; n1 < N; n1++)				
			for(int n2 = 0; n2 < N; n2++)						
				if(n1 != n2)				
					traffMatrixTelcoTelco.set(n1, n2, popuWeightVector[n1]*popuWeightVector[n2]);					
		
		DoubleMatrix2D initialTraffMatrixTelcoTelco = DoubleFactory2D.dense.make(N,N,0);							
		DoubleMatrix2D traffMatrix = DoubleFactory2D.dense.make(N,N,0);		
		List<DoubleMatrix2D> traffMatrixPerCDN = new ArrayList<DoubleMatrix2D>();
		
		double[] numberOfAccesses = appAndCDNInfo.computeNumberOfAccesses(U);
		List<List<List<List<Integer>>>> replicaPlacements = appAndCDNInfo.computeInitialReplicaPlacements(cdnIndexesNodes);			
		
		for (int y = 0; y < simYears; y++)
		{
			
//			double timenow = (System.nanoTime()-iniTime)*1e-9; 
//			System.out.println("Year "+y+ " in " + timenow + " seconds.");
			double[] carriedTrafficD2CPerService = new double[S];
			double[] carriedTrafficD2DPerService = new double[S];
			double[] offeredTrafficD2CPerService = new double[S];
			int[] numberOfCarriedDemandsPerService = new int[S];
			double[] propagationTime = new double[S];
			double[] rttPerService = new double[S];				

			netPlan = originalnetPlan.copy();
			List<Node> netNodes  = netPlan.getNodes();
			List<Link> netLinks = netPlan.getLinks();			
			
			// Telco-Telco Traffic Matrix
			telecoTelcoTraffic = H_TelcoTelco*Math.pow((1+CAGR_telcoTelco), y);
			initialTraffMatrixTelcoTelco = TrafficMatrixGenerationModels.normalizationPattern_totalTraffic(traffMatrixTelcoTelco, telecoTelcoTraffic);
			List<Demand> notCDNDmands = netPlan.addDemandsFromTrafficMatrix(initialTraffMatrixTelcoTelco);			
			
			for (Demand d: notCDNDmands)
			{
				seqLinks = GraphUtils.getShortestPath(netNodes,netLinks, d.getIngressNode(), d.getEgressNode(), null);
				netPlan.addRoute(d, d.getOfferedTraffic(), d.getOfferedTraffic(), seqLinks, null);
			}

			stat_carriedNotCDNTraffic.add(y,netPlan.getDemandTotalCarriedTraffic());
			
			// Compute traffic per application
			traffMatrixPerCDN.clear();
			for (int c = 0; c < C; c++)
				traffMatrixPerCDN.add(c, DoubleFactory2D.sparse.make(N,N));
			
//			netPlan.removeAllDemands();
			netPlan = originalnetPlan.copy();
			netNodes  = netPlan.getNodes();
			netLinks = netPlan.getLinks();
			
			int appIndex = 0;
			for (Triple<Double,Integer, int[]> app : appInfo)
			{						
				appService = app.getSecond();														// Type of service of the app
				cagr_service = servicesInfo.get(appService).getSecond();
				beta_s = servicesInfo.get(appService).getThird();
				appTraffic = H_D2C*app.getFirst()*Math.pow((1+cagr_service),y);
				int[] appCDNs = app.getThird();												      	// List of CDNs where the app can carry the traffic				
				numberOfCDNs = appCDNs.length;																	// Number of different CDN available in this app
				h_ac = appTraffic/(double) numberOfCDNs;											// Traffic per each CDN for the app					
				
//				List<Integer> cdnIndexNodes = cdnIndexesNodes.get(appService);
				Set<Node> destinationNodes = new HashSet<Node>();
				DoubleMatrix2D traffMatrixAppCDN = DoubleFactory2D.dense.make(N,N);			
				
				for (int c=0; c < numberOfCDNs; c ++)
				{
					for (int i=0; i < U; i++)					
					{
						// D2C Traffic Matrices
	
	//					Nc = uPlacementIndexes.size();
						List<Integer> uPlacementIndexes = replicaPlacements.get(appIndex).get(c).get(i);
						traffMatrixAppCDN = DoubleFactory2D.dense.make(N,N);
						
						for(int n1 = 0; n1 < N; n1++)				
							for(int n2 = 0; n2 < N; n2++)						
								if(n1 != n2)
									if(uPlacementIndexes.contains(n1))								
										traffMatrixAppCDN.set(n1, n2, popuWeightVector[n1]*popuWeightVector[n2]);										
						
						traffMatrix = TrafficMatrixGenerationModels.normalizationPattern_totalTraffic(traffMatrixAppCDN, h_ac*numberOfAccesses[i]);		
						traffMatrixPerCDN.get(appCDNs[c]).assign(traffMatrix,DoubleFunctions.plus);
						
						for (Demand d : netPlan.addDemandsFromTrafficMatrix(traffMatrix))
						{								
							offeredTrafficD2CPerService[appService] += d.getOfferedTraffic();
							numberOfCarriedDemandsPerService[appService] ++;
							if (rand.nextDouble() <= 0.8) 	
							{						
								if (!uPlacementIndexes.contains(d.getEgressNode().getIndex()) )						
								{
									seqLinks = GraphUtils.getShortestPath(netNodes, netLinks, d.getIngressNode(), d.getEgressNode(), null);
									Route r = netPlan.addRoute(d, d.getOfferedTraffic(), d.getOfferedTraffic(), seqLinks, null);
									carriedTrafficD2CPerService[appService] += d.getOfferedTraffic();
									propagationTime[appService] += r.getPropagationDelayInMiliseconds();
								}							
							}
							else
							{
								seqLinks  = GraphUtils.getShortestPath(netNodes, netLinks, d.getIngressNode(), d.getEgressNode(), null);
								Route r = netPlan.addRoute(d, d.getOfferedTraffic(), d.getOfferedTraffic(), seqLinks, null);
								carriedTrafficD2CPerService[appService] += d.getOfferedTraffic();
								propagationTime[appService] += r.getPropagationDelayInMiliseconds();
							}				
						}
					
						// D2D Traffic Matrices (InterCDN Traffic)		
					
						Nmas = numberOfAccesses[i];
						int Nr = uPlacementIndexes.size();
						replicaTraffic = beta_s*h_ac*Nmas;	
	
						destinationNodes.clear();
						originNodeReplicaIndex = uPlacementIndexes.get(0);
						Node originNode = netPlan.getNode(originNodeReplicaIndex);
						
						for(int r = 1; r < Nr; r ++)
							if (uPlacementIndexes.get(r) != originNodeReplicaIndex)
								destinationNodes.add(netPlan.getNode(uPlacementIndexes.get(r)));	
						
						if(!destinationNodes.isEmpty())
						{
							multicastDemand = netPlan.addMulticastDemand(originNode, destinationNodes, replicaTraffic, null);						
							multLinks = GraphUtils.getMinimumCostMulticastTree(netLinks, netPlan.getMatrixNodeLinkOutgoingIncidence(), netPlan.getMatrixNodeLinkIncomingIncidence(),
									linkCost, originNode, destinationNodes, -1, -1, -1, -1, "cplex", path , 10);						
							multicastTree = netPlan.addMulticastTree(multicastDemand, replicaTraffic, replicaTraffic, multLinks, null);
							carriedTrafficD2DPerService[appService] += multicastTree.getCarriedTraffic();
						}
					}		
				}
				appIndex++;
			}
					
			for(int s=0; s<S; s++)			
				rttPerService[s] = propagationTime[s]/(double) numberOfCarriedDemandsPerService[s];			
			
			stat_carried.add(y, netPlan.getDemandTotalCarriedTraffic());
			stat_offered.add(y, netPlan.getDemandTotalOfferedTraffic());
			stat_multiCastTraffic.add(y, netPlan.getMulticastDemandTotalCarriedTraffic());
			stat_avRTT.add(y, rttPerService);
			stat_carriedD2CTrafficPerService.add(y, carriedTrafficD2CPerService);	
			stat_offeredD2CTrafficPerService.add(y, offeredTrafficD2CPerService);
			int newDC = 0;	
			// Update Data Center Locations
			if(G.getDouble() > 0)	
				for (int a = 0; a < A; a++)
					for(int c = 0; c < appInfo.get(a).getThird().length; c++)
					{
						int[] thisCDN = appInfo.get(a).getThird();
						newDC += appAndCDNInfo.addDcIntoCDN(netPlan, a, thisCDN[c], G.getDouble(),traffMatrixPerCDN.get(thisCDN[c]));
					}
			
			stat_numberOfNewDC.add(y+1,newDC);
			replicaPlacements = appAndCDNInfo.checkReplicaPlacements();			
		}
					
		String root;
		if (isLocal.getBoolean()) root = "C:/Users/Javi/OneDrive/Projects/Proyecto Trend CDN/Results/";
		else root = "../trendTraffic/Results/";
		
		String gString = Double.toString(G.getDouble());
		String simString = Integer.toString(sim.getInt());
		
		stat_carried.printToFile(new File(root+"carried"+gString+"_"+simString+".txt"));
		stat_offered.printToFile(new File(root+"offered"+gString+"_"+simString+".txt"));
		stat_multiCastTraffic.printToFile(new File(root+"multicast"+gString+"_"+simString+".txt"));
		stat_avRTT.printToFile(new File(root+"avRtt_"+gString+"_"+simString+".txt"));		
		stat_carriedNotCDNTraffic.printToFile(new File(root+"carriedNotCDNTraffic_"+gString+"_"+simString+".txt"));		
		stat_carriedD2CTrafficPerService.printToFile(new File(root+"carriedD2CTrafficPerService_"+gString+"_"+simString+".txt"));
		stat_offeredD2CTrafficPerService.printToFile(new File(root+"offeredD2CTrafficPerService_"+gString+"_"+simString+".txt"));
		stat_numberOfNewDC.printToFile(new File(root+"numberOfNewDC"+gString+"_"+simString+".txt"));
		
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
	
}
