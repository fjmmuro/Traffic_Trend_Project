cd c:/Net2Plan-0.4.2

java -jar Net2Plan-cli.jar --mode net-design ^
--input-file C:\net2plan-0.4.2\workspace\data\networkTopologies\NSFNet_N14_E42_highcapacity.n2p ^
--class-file C:\Users\Javi\Documents\Traffic_Trend_Project\bin\com\net2plan\general\trafficTrendCDN.class ^
--class-name trafficTrendCDN ^
--output-file C:\Users\Javi\OneDrive\Projects\output.n2p ^
--alg-param A=10 ^
--alg-param sim=5 ^
--alg-param simYears=5 ^
--alg-param isLocal=true

