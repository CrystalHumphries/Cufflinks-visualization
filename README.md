Cufflinks-Visualization
=======================

Purpose: created BED files to visualize different transcripts on USSC browser from cufflinks output

How to Run:
  perl BEDfromTrancripts.pl -T transcripts_file.txt -L gencode.attributes.file.txt

The transcripts file is a typical output from cuffcompare:

8\_00458540	XLOC\_011348	FBXO2|FBXO2	x	-	-	q3:AD2.974|AD2.974.1|100|6.152141|4.722721|7.581561|58.607883|430	q4:AD3.95|AD3.95.1|100|9.046277|5.273719|12.818835|21.328327|395	q5:AD4.565|AD4.565.1|100|9.136013|7.403165|10.868861|78.222360|461	q6:AD5.1346|AD5.1346.1|100|11.989985|10.757990|13.221980|128.213270|672	q7:AD6.286|AD6.286.1|100|9.181492|6.215037|12.147948|24.486141|438	q8:AD7.429|AD7.429.1|100|10.805550|9.351270|12.259830|77.007361|696	q9:AD8.809|AD8.809.1|100|13.744011|11.935960|15.552061|103.957922|568	q10:AD9.446|AD9.446.1|100|10.275882|8.469653|12.082110|51.629598|676	q11:AD10.329|AD10.329.1|100|2.582497|1.617748|3.547246|18.637065|360	q12:CT1.436|CT1.436.1|100|11.942306|10.506601|13.378011|125.901033|609	q13:CT2.840|CT2.840.1|100|10.398910|8.814966|11.982853|166.329675|456	q14:CT3.479|CT3.479.1|100|4.004266|3.202365|4.806167|59.372393|435	q15:CT4.346|CT4.346.1|100|6.608667|5.253728|7.963606|69.597980|365	q16:CT5.573|CT5.573.1|100|12.035976|10.546960|13.524992|110.168747|680	-	q18:CT7.491|CT7.491.1|100|16.836556|13.956757|19.716356|110.514101|358	q19:CT8.464|CT8.464.1|100|36.136830|29.070582|43.203078|172.204315|299	q20:CT9.459|CT9.459.1|100|11.326804|8.865338|13.788270|84.022312|438	-	q22:D1.447|D1.447.1|100|7.932166|6.075806|9.788525|60.498912|423	q23:D2.1805|D2.1805.1|100|2.149630|1.586002|2.713258|28.950735|449	q24:D3.1087|D3.1087.1|100|44.130410|38.356418|49.904403|277.076375|307	-	q26:D5.610|D5.610.1|100|7.096536|5.524014|8.669059|49.536377|417	-	q28:D7.730|D7.730.1|100|7.926554|6.696342|9.156766|59.980629|630	q29:D8.315|D8.315.2|55|1.928300|1.209276|2.647323|13.749436|629	q30:D9.482|D9.482.1|100|3.177128|2.382846|3.971410|20.911024|572	q31:D10.233|D10.233.1|100|7.543618|5.430982|9.656254|18.714736|572

The output is in a BED file which denotes exon-intron junctions.
In addition, the transcript is color coded based on which group it is expressed in. 

	  #lime   if present in CT only 
	  #blue   if present in AD + CT only
	  #yellow if present in DLB only
	  #purple if present in AD + DLB 
	  #teal   if present in DLB + CT
	  #black  if present in all three groups

Example Output
-------------------
chr1	11637558	11644939	8\_00433015\_AD1\_CT2	500	+	11637558	11644939	0,0,255	5	38,294,97,135,1165	0,892,3353,3822,6216

How to Visualize with UCSC Genome Browser
-----------------------------------------------
	1. go to this website <http://genome.ucsc.edu/cgi-bin/hgCustom?hgHubConnect.destUrl=..%2Fcgi-bin%2FhgTracks&clade=mammal&org=Human&db=hg18&position=chr19%3A996%2C074-996%2C679&hgt.positionInput=enter+position%2C+gene+symbol+or+search+terms&hgt.suggestTrack=knownGene&hgsid=363174797> 
	2. place output
	3. click submit
	4. click on hyperlink under "pos"
	5. zoom out to view full transcript
