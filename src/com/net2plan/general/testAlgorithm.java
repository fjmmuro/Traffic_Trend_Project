package com.net2plan.general;


import cern.colt.matrix.tdouble.DoubleFactory1D;
import cern.colt.matrix.tdouble.DoubleFactory2D;
import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.jet.math.tdouble.DoubleFunctions;
import com.net2plan.interfaces.networkDesign.*;
import com.net2plan.libraries.GraphUtils;
import com.net2plan.libraries.TrafficMatrixGenerationModels;
import com.net2plan.utils.Constants.RoutingType;
import com.net2plan.utils.*;

import java.io.File;
import java.util.*;

/** This is a template to be used in the lab work, a starting point for the students to develop their programs
 * 
 */
public class testAlgorithm implements IAlgorithm
{

    private InputParameter simYears = new InputParameter ("simYears", (int) 20 , "Simulation Years from 2015");
    private InputParameter G = new InputParameter ("G", (double) 1 , "Expansion Factor per CDN");
    private InputParameter sim = new InputParameter ("sim", (int) 1, "simulation index");
    private InputParameter isLocal = new InputParameter("isLocal", (boolean) false, "false if the execution is run in a server, true otherwise");
    private InputParameter A = new InputParameter("A", (int) 5, "Number of Apps");
    private InputParameter S = new InputParameter("S", (int) 5, "Number of Services");
    private InputParameter C = new InputParameter("C", (int) 5, "Number of CDNs");
    private InputParameter U = new InputParameter("U", (int) 100, "Number of Content Units");
    private InputParameter avNumReplicasPerContentUnit = new InputParameter("avNumReplicasPerContentUnit", (double) 2, "Average number of replicas of a content unit");
    private InputParameter rngSeed = new InputParameter("rngSeed", (long) 0, "The seed for the random numbers generator");
    private InputParameter debugMode = new InputParameter("debugMode", (boolean) true, "Set up true if debug mode");

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
	
	@Override
	public String executeAlgorithm(NetPlan netPlan, Map<String, String> algorithmParameters, Map<String, String> net2planParameters)
	{
        long iniTime = System.nanoTime();
        netPlan.removeAllDemands();

		/* Initialize all InputParameter objects defined in this object (this uses Java reflection) */
        InputParameter.initializeAllInputParameterFieldsOfObject(this, algorithmParameters);

        // Set the population of the nodes

        // Initial checks
        final int N = netPlan.getNumberOfNodes();
        final int E = netPlan.getNumberOfLinks();
        final int S = this.S.getInt();
        final int C = this.C.getInt();
        final int A = this.A.getInt();
        final int U = this.U.getInt();
        final int simYears = this.simYears.getInt();
        final double H_D2C = 0.76*1000;
        final double H_TelcoTelco = 1000-H_D2C;
        final Random rand = new Random(this.rngSeed.getLong ());

        if (netPlan.getNode(0).getAttribute("population") == null)
            TrafficTrendUtils.setNodePopulation(netPlan);

        double[] populationWeightVector = new double[N];

        final double CAGR_telcoTelco = 0.1;
        netPlan.setRoutingType(RoutingType.SOURCE_ROUTING);

        // Read services, applications and available CDNs
        TrafficTrendUtils appAndCDNInfo = new TrafficTrendUtils(C,S,A,U,rand);
        appAndCDNInfo.generateServicePerApp();

        String path = "";
        if (!isLocal.getBoolean())
            path = "/opt/ibm/ILOG/CPLEX_Studio127/cplex/bin/x86-64_linux/libcplex1270.so";

        List<Triple<Double, Integer, int[]>> appInfo = appAndCDNInfo.getApps();
        List<List<Node>> nodesWithDCPerCDN_c = appAndCDNInfo.getInitialDCPlacementPerCDN(netPlan , C , rand);

        int cdn = 0;
        for (List<Node> nodesThisCDN : nodesWithDCPerCDN_c)
        {
            System.out.println(" ---- CDN " + cdn + " ----");
            for (Node node : nodesThisCDN)
                System.out.print(" " + node.getIndex());
            cdn++;
            System.out.println(" ");
        }

        // Services Test

		List<Triple<Double,Double,Double>> servicesInfo = appAndCDNInfo.getServices();

        int serviceID = 0;
        for ( Triple<Double,Double,Double> service : servicesInfo)
        {
            System.out.println(" ---- Service " + serviceID + " -----");
            System.out.println(" * Traffic: " + service.getFirst());
            System.out.println(" * CAGR: " + service.getSecond());
            System.out.println(" * Beta: " + service.getThird());
            serviceID++;
        }

        // Initial Replica Placements Test

        final DoubleMatrix1D populationPerNode = NetPlan.getAttributeValues(netPlan.getNodes() , "population" , 0);
        final double totalPopulation = populationPerNode.zSum();
        for (int n = 0; n< N; n++)
            populationWeightVector[n] = populationPerNode.get(n)/totalPopulation;

        final Map<Pair<Node,Node>,List<List<Link>>> cpl = netPlan.computeUnicastCandidatePathList (null , 1, -1, -1, -1, -1,-1, -1, null);
		final DoubleMatrix2D rtt_n1n2 = DoubleFactory2D.dense.make(N,N);
		final DoubleMatrix2D numHops_n1n2 = DoubleFactory2D.dense.make(N,N);
		for (Pair<Node,Node> n1n2 : cpl.keySet())
		{
			final List<Link> seqLinks = cpl.get(n1n2).get(0);
			final double rtt = seqLinks.stream().mapToDouble(e->e.getPropagationDelayInMs()).sum();
			rtt_n1n2.set(n1n2.getFirst().getIndex() , n1n2.getSecond().getIndex() , rtt);
			numHops_n1n2.set(n1n2.getFirst().getIndex() , n1n2.getSecond().getIndex() , seqLinks.size());
		}
		List<List<List<Node>>> replicaPlacements = appAndCDNInfo.computeReplicaPlacementsForAllCDNs (netPlan , avNumReplicasPerContentUnit.getDouble() , numHops_n1n2 , populationWeightVector, path, true);

//		int app = 0;
//		for(List<List<Node>> replicasThisCDN : replicaPlacements)
//		{
//		    int thisCdn = 0;
//            System.out.println(" ");
//            System.out.println(" ---- CDN " + thisCdn + " ----");
//            for (List<Node> replicaThisAppThisCDN : replicasThisCDN)
//            {
//                int cu = 0;
//                for (Node replica : replicaThisAppThisCDN) System.out.print(" " + replica.getIndex());
//                cu++;
//
//                System.out.println(" ");
//                thisCdn++;
//            }
//            app++;
//        }


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
		final List<Demand> notCDNDemands = netPlan.addDemandsFromTrafficMatrix(initialTraffMatrixTelcoTelco);

		for (Demand d: notCDNDemands)
		{
			final List<Link> seqLinks = GraphUtils.getShortestPath(netPlan.getNodes() ,netPlan.getLinks(), d.getIngressNode(), d.getEgressNode(), null);
			netPlan.addRoute(d, d.getOfferedTraffic(), d.getOfferedTraffic(), seqLinks, null);
		}
		return Pair.of(netPlan.getVectorDemandOfferedTraffic().zSum() , netPlan.getVectorLinkCarriedTraffic().zSum());
	}
	
}
