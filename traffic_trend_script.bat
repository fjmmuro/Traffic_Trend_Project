cd c:/net2plan-0.5-SNAPSHOT

for %%g in (0 1 2) do (
	for %%s in (1 2 3 4 5) do (
		java -jar Net2Plan-cli.jar --mode net-design ^
		--input-file C:\net2plan-0.5-SNAPSHOT\workspace\data\networkTopologies\NSFNet_N14_E42.n2p ^
		--class-file C:\Users\javie\Documents\Git\Traffic_Trend_Project\out\production\Traffic_Trend_Project\com\net2plan\general\TrafficTrendCDN.class ^
		--class-name TrafficTrendCDN ^
		--output-file C:\Users\javie\OneDrive\Projects\output.n2p ^
		--alg-param A=5 ^
		--alg-param U=100 ^
		--alg-param G=%%g ^
		--alg-param rngSeed=0 ^
		--alg-param sim=%%s ^
		--alg-param avNumReplicasPerContentUnit=2.0 ^
		--alg-param simYears=15 ^
		--alg-param isLocal=true
	)

)
