import java.io.File;
import java.io.FileNotFoundException;
import java.util.Arrays;
import java.util.Random;
import java.util.Scanner;

class pBlast {
	public static void main(String[] args) {
		// Read FASTA file
		String fName = "E:\\Sylvia\\Workspace\\FinalProject\\src\\dna.fa";
		String dna="";
		try {
			Scanner sc = new Scanner(new File(fName));
			while (sc.hasNextLine()) {
				dna = dna + sc.nextLine();				
			}
			sc.close();
		}
		catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		
		fName = "E:\\Sylvia\\Workspace\\FinalProject\\src\\probability";
		double[] prob = new double[dna.length()];
		try {
			Scanner sc = new Scanner(new File(fName));
			int count = 0;
			while (sc.hasNextDouble()) {
				prob[count] = sc.nextDouble();
				count++;
			}
			sc.close();
		}
		catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		
		
		// Set query
		String query = randomize(generate(dna, prob, 1000, 5000), 0.05); 
		
		
		String[][] preData = preProcess(dna, prob);
		
		String[][] seed = findSeed(query, preData);
		
		String[] hsp = ungap(seed, query, dna, prob);
		
		// find the best hsp
		String[] index = hsp[0].split("[;]+");
		double[] score = convScore(hsp[1].split("[ ]+"));
		
		
		double threshold = 0.1 * query.length();
		String[] bestAlign = new String[3];
		bestAlign[2] = "0";
		for(int i = 0; i< score.length;i++) {
			if (score[i] > threshold) {
				String[] ind = index[i].split("[ ]+");
				String[] nAlign = gapped(ind[0],ind[1],ind[2],ind[3], query, dna, prob, score[i]);
				if (Double.parseDouble(bestAlign[2]) < Double.parseDouble(nAlign[2])) {
					bestAlign= nAlign;
					//[0] = nAlign[0];
					//bestAlign[1] = nAlign[1];
					//bestAlign[2] = nAlign[2];
				}
			}
		}
		
		
		System.out.println(query);
		System.out.println(bestAlign[0]);
		System.out.println(bestAlign[1]);
		System.out.println(bestAlign[2]);
		System.out.println(bestAlign[3]);
		System.out.println(bestAlign[4]);
		
	}
	
	public static String randomize(String q, double percent) {
		String query = "";
		int len = q.length();
		Random index = new Random();
		String nucleotide = "ACGT";
		int[] indel = new int[(int) Math.floor(len*percent)];
		for(int i = 0; i < indel.length; i++) {
			indel[i] = index.nextInt(len);
		}
		
		Arrays.sort(indel);
		int prev = 0;
		for(int i = 0; i< indel.length; i++) {
			double r = Math.random();
			if (r <= 1/3) { // deletion
				if(prev >= indel[i]) {
					prev = indel[i] + 1;
					continue;
				}
				else{
					query = query + q.substring(prev, indel[i]);
					prev = indel[i] + 1;
				}
				
			}
			else if (r <= 2/3) { //insertion
				if(prev >= indel[i]) {
					prev = indel[i];
					continue;
				}
				else {
					query  = query + q.substring(prev, indel[i]);
					prev = indel[i];
					query = query + nucleotide.charAt(index.nextInt(4));
				}
				
			}
			else {
				if(prev >= indel[i]) {
					prev = indel[i]+1;
					continue;
				}
				else {
					query = query + q.substring(prev, indel[i]);
					prev = indel[i] + 1;
					query = query + nucleotide.charAt(index.nextInt(4));
				}
			}
		}
		if(prev < len) {
			query = query + q.substring(prev, len);
		}
		return query;
	}
	
	
	
	public static String generate(String dna, double[] prob, int pos, int len) {
		String query = "";
		for (int i = 0; i<len; i++) {
			double r = Math.random();
			double p = prob[pos+i];
			String c = "ACGT";
			if(c.indexOf(dna.charAt(pos+i))==0) {
				c = "CGT";
			}
			else if(c.indexOf(dna.charAt(pos+i))==1){
				c = "AGT";
			}
			else if(c.indexOf(dna.charAt(pos+i))==2){
				c = "ACT";
			}
			else{
				c = "ACG";
			}
			
			if(r <= p) {
				query = query + dna.charAt(pos+i);
			}
			else if(r <= (p+(1-p)/3)) {
				query = query + c.charAt(0);
			}
			else if(r <= (p + 2*(1-p)/3)) {
				query = query + c.charAt(1);
			}
			else {
				query = query + c.charAt(2);
			}
		}
		
		return query;
	}
	
	public static String[] gapped(String qS, String qE, String dS, String dE, String query, String dna, double[] prob, double score) {
		String[] align = new String[5];
		
		int qStart = Integer.parseInt(qS);
		int dStart = Integer.parseInt(dS);
		int qEnd = Integer.parseInt(qE);
		int dEnd = Integer.parseInt(dE);
		
		int dRange;
		if(dStart+1 > 4*(qStart+1)) {
			dRange = 4*(qStart+1);
		}
		else {
			dRange = dStart+1 ;
		}
		
		double[][] lMatch = new double[qStart+1][dRange];
		int[][] lPointer = new int[qStart+1][dRange];
		String lq = "";
		String ld = "";
		double[] best = {0,0,0};
		
		// fill in left matrix
		for(int i=0;i<qStart+1;i++) {
			for(int j = 0; j<dRange;j++) {
				//indel penalty = -1
				if(i==0 && j==0) {
					lMatch[i][j] = 0;
				}
				else if(i==0) {
					lMatch[i][j] = lMatch[i][j-1] -1;
					lPointer[i][j] = 2;
				}
				else if(j==0) {
					lMatch[i][j] = lMatch[i-1][j] -1;
					lPointer[i][j] = 1;
				}
				else {
					double m = lMatch[i-1][j-1] + match(query.charAt(qStart-i), dna.charAt(dStart-j), prob[dStart-j]);
					double iDel = lMatch[i-1][j] -1;
					double jDel = lMatch[i][j-1] -1;
					
					if(m>iDel && m> jDel) {
						lMatch[i][j] = m;
						lPointer[i][j] = 0;
					}
					else if(iDel > m && iDel > jDel) {
						lMatch[i][j] = iDel;
						lPointer[i][j] = 1;
					}
					else {
						lMatch[i][j] =jDel;
						lPointer[i][j] = 2;
					}
					
					if(lMatch[i][j] > best[0]) {
						best[0] = lMatch[i][j];
						best[1] = i;
						best[2] = j;
					}
				}
			}
		}	
		
		
		int x = (int)best[1];
		int y = (int)best[2];
		align[3] = (qStart - x) +" ";
		align[4] = (dStart - y) +" ";
		// Trace back
		while (x>0 || y>0) {
			if(lPointer[x][y] == 0) {
					lq = lq + query.charAt(qStart-x);
					ld = ld + dna.charAt(dStart-y);
					x--;
					y--;
				}
			else if(lPointer[x][y] == 1) {
					lq = lq + query.charAt(qStart-x);
					ld = ld + "-";
					x--;
			}
			else {
					lq = lq + "-";
					ld = ld + dna.charAt(dStart-y);
					y--;
			}
		}
		score = score + best[0];
		
		
		
		// fill in right matrix
		best[0] = 0;
		best[1] = 0;
		best[2] = 0;
		
		if(dna.length()-dEnd > 4*(query.length()-qEnd)) {
			dRange = 4*(query.length()-qEnd);
		}
		else {
			dRange = dna.length()-dEnd;
		}
		
		double[][] rMatch = new double[query.length()-qEnd][dRange];
		int [][] rPointer = new int[query.length()-qEnd][dRange];
		String rq = "";
		String rd = "";
		
		
		for(int i=0;i<query.length()-qEnd;i++) {
			for(int j = 0; j<dRange;j++) {
				//indel penalty = -1
				if(i==0 && j==0) {
					rMatch[i][j] = 0;
				}
				else if(i==0) {
					rMatch[i][j] = rMatch[i][j-1] -1;
					rPointer[i][j] = 2;
				}
				else if(j==0) {
					rMatch[i][j] = rMatch[i-1][j] -1;
					rPointer[i][j] = 1;
				}
				else {
					double m = rMatch[i-1][j-1] + match(query.charAt(qEnd+i), dna.charAt(dEnd+j), prob[dEnd+j]);
					double iDel = rMatch[i-1][j] -1;
					double jDel = rMatch[i][j-1] -1;
					
					if(m>iDel && m> jDel) {
						rMatch[i][j] = m;
						rPointer[i][j] = 0;
					}
					else if(iDel > m && iDel > jDel) {
						rMatch[i][j] = iDel;
						rPointer[i][j] = 1;
					}
					else {
						rMatch[i][j] =jDel;
						rPointer[i][j] = 2;
					}
					
					if(rMatch[i][j] > best[0]) {
						best[0] = rMatch[i][j];
						best[1] = i;
						best[2] = j;
					}
				}
			}
		}	
		
		// Trace back
		x = (int)best[1];
		y = (int)best[2];
		align[3] = align[3] + (qEnd + x);
		align[4] = align[4] + (dEnd + y);
		
		while(x > 0 || y > 0) {
			if(rPointer[x][y] == 0) {
				rq = query.charAt(qEnd+x) +rq;
				rd = dna.charAt(dEnd+y) +rd;
				x--;
				y--;
			}
			else if(rPointer[x][y] == 1) {
				rq = query.charAt(qStart+y) +rq;
				rd = "-" +rd;
				x--;
			}
			else {
				rq = "-" +rq;
				rd = dna.charAt(dStart+y) +rd;
				y--;
			}	
		}
		score = score + best[0];
		
		
		align[0] = lq + query.substring(qStart, qEnd+1) + rq;
		align[1] = ld + dna.substring(dStart, dEnd+1) + rd;
		align[2] = "" + score;
	
		return align;
	}
	
	private static double[] convScore(String[] sStr) {
		double[] score = new double[sStr.length];
		
		for(int i = 0;i<sStr.length;i++) {
			score[i] = Double.parseDouble(sStr[i]);
		}
		
		return score;
	}
	
	
	public static String[] ungap(String[][] seed, String query, String dna, double[] prob) {
		//String[][] hsp = new String[seed.length][3];
		String[] hsp = {"",""};
		
		int qStart, qEnd, dStart, dEnd, left, right, maxL, maxR;
		double score, max, t;
		
		t = 5.0; // threshold
		
		for(int i = 0; i<seed.length; i++) {
			if(seed[i][0] == null) {
				continue;
			}
			
			String[] dPos = seed[i][0].split("[ ]+");
			String[] s = seed[i][1].split("[ ]+");
			
			for(int j=0; j<dPos.length;j++) {
				qStart = i;
				qEnd = i+10;
				
				dStart = Integer.parseInt(dPos[j]);
				dEnd = dStart + 10;
				
				score = Double.parseDouble(s[j]);
				max = score;
				
				left = 1;
				maxL = 0;
				// Extend to left
				while((qStart-left)>=0 && (dStart-left)>=0) {
					score = score + match(query.charAt(qStart-left), dna.charAt(dStart-left), prob[dStart-left]);
					if(score >= max) {
						max = score;
						maxL = left;
					}
					left++;
					if((max-score) > t) {
						break;
					}
				}
				
				score = max;
				qStart = qStart - maxL;
				dStart = dStart - maxL;
				
				right = 1;
				maxR = 0;
				// Extend to right
				while((qEnd+right) < query.length() && (dEnd+right) < dna.length()) {
					score = score + match(query.charAt(qEnd+right), dna.charAt(dEnd+right), prob[dEnd+right]);
					
					if(score >= max) {
						max = score;
						maxR = right;
					}
					right++;
					if((max-score) > t) {
						break;
					}
				}
				
				score = max;
				qEnd = qEnd + maxR;
				dEnd = dEnd + maxR;
				
				String nHsp = qStart + " " + qEnd + " " + dStart + " " + dEnd + ";";
				
				// avoid repeating
				if(hsp[0].indexOf(nHsp) < 0) {
					hsp[0] = hsp[0] + nHsp;
					hsp[1] = hsp[1] + score + " ";
				}
			}
		}
		
		
		return hsp;
	}
	
	private static double match(char q, char d, double p) {
		if(q==d) {
			return p;
		}
		else {
			return -(p-(1-p)/3);
		}
	}
	
	public static String[][] findSeed(String q, String[][] preData) {
		String[][] seed = new String[q.length()-10][2];
		
		int[] query = convert(q);
		
		int qPos = 0;
		for(int i=0; i<11; i++) {
			qPos = qPos*4 + query[i];
		}
		
		seed[0][0] = preData[qPos][0];
		seed[0][1] = preData[qPos][1];
		
		int fir = query[0];
		
		for(int i=11; i<query.length;i++) {
			if(fir != 0) {
				qPos = (qPos - fir*(1048576))*4 + query[i];
				seed[i-10][0] = preData[qPos][0];
				seed[i-10][1] = preData[qPos][1];
			}
			else {
				qPos = qPos*4 + query[i];
				seed[i-10][0] = preData[qPos][0];
				seed[i-10][1] = preData[qPos][1];
			}
			fir = query[i-10];
		}
		
		
		return seed;
	}
	
	
	public static int[] convert(String o) {
		int[] n = new int[o.length()];
		
		for(int i=0; i<o.length();i++) {
			if(o.charAt(i) == 'A') {
				n[i] = 0;
			}
			else if(o.charAt(i)  == 'C') {
				n[i] = 1;
			}
			else if(o.charAt(i)  == 'G') {
				n[i] = 2;
			}
			else {
				n[i] = 3;
			}
		}
		
		return n;
	}
	
	
	public static String[][] preProcess(String dna, double[] prob){
		String[][] pre = new String[4194304][2];
		
		//set seed threshold
		double t = 9;
		
		
		//convert to base 4 digit
		int[] dnaNum = convert(dna);
		
		//find the first seed position
		int pos = 0;
		double conf = 0;
		int fir = dnaNum[0];
		double prev = prob[0];
		
		for(int i=0;i<11;i++) {
			pos = pos*4+dnaNum[i];
			conf = conf  + prob[i];
		}
		
		if (conf > t) {
			pre[pos][0] = "0 ";
			pre[pos][1] = conf + " ";
		}
		
		//find all seed position
		for(int i = 11; i< dnaNum.length;i++) {
			conf = conf - prev + prob[i];
			prev = prob[i-10];
			
			if(fir != 0) {
				pos = (pos - fir*(1048576))*4 + dnaNum[i];
			}
			else {
				pos = pos*4 + dnaNum[i];	
			}
			
			if (conf > t) {
				if (pre[pos][1]== null) {
					pre[pos][1] = conf + " ";
				}
				else {
					pre[pos][1] = pre[pos][1] + conf + " ";
				}
				
				if (pre[pos][0]==null) {
					pre[pos][0] = (i-10) + " ";
					
				}
				else {
					pre[pos][0] = pre[pos][0] + (i-10) + " ";
				}
			}
			
			fir = dnaNum[i-10];
			
		}
				
		return pre;
	}
}
