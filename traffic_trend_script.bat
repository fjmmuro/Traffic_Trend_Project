cd c:/net2plan-0.5-SNAPSHOT

java -jar Net2Plan-cli.jar --mode net-design ^
--input-file C:\net2plan-0.5-SNAPSHOT\workspace\data\networkTopologies\NSFNet_N14_E42.n2p ^
--class-file C:\Users\Javi\Documents\Traffic_Trend_Project\bin\com\net2plan\general\TrafficTrendCDN.class ^
--class-name TrafficTrendCDN ^
--output-file C:\Users\Javi\OneDrive\Projects\output.n2p ^
--alg-param A=10 ^
--alg-param U=100 ^
--alg-param rngSeed=0 ^
--alg-param sim=2 ^
--alg-param avNumReplicasPerContentUnit=2.0 ^
--alg-param simYears=20 ^
--alg-param isLocal=true

