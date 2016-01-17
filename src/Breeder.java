import java.util.*;

public class Breeder{
	
	private int num;
	private int speciesNum;
	private int nodes;
	private int innovTotal;
	private double delta = 10;
	private HashMap<Integer, Gene> geneList;
	private HashMap<Integer, HashMap<Integer, Gene>> generation;
	private HashMap<Integer, Integer> species; //The key is the organism number, the value is the species number
	private HashMap<Integer, Double> fitness; //The key is organism number, the value is fitness number
	private Game game;
	
	public Breeder(int num){
		
		this.num = num;
		speciesNum = 0;
		nodes = 20;
		innovTotal = 0;
		
		geneList = new HashMap<Integer, Gene>();
		generation = new HashMap<Integer, HashMap<Integer, Gene>>();
		species = new HashMap<Integer, Integer>();
		fitness = new HashMap<Integer, Double>();
		
		game = new Game();
		
		makeFirstGen();
		
	}
	
	private void makeFirstGen(){
		
		//Create all possible starting genes
		for(int i = 0; i < 4; i++){
			
			for(int j = 0; j < 16; j++){
				
				//even genes are positive connections
				Gene gene = new Gene(i, j + 4, 1.0, innovTotal);
				geneList.put(innovTotal, gene);
				innovTotal++;
				//odd genes are negative
				Gene gene2 = new Gene(i, j + 4, -1.0, innovTotal);
				geneList.put(innovTotal, gene2);
				innovTotal++;
				//let's stick to even/odd
				
			}
			
		}
		
		Random random = new Random();
		
		for(int i = 0; i < num; i++){
			
			int numGenes = 32 + (int)(16.0*random.nextGaussian());
			
			if(numGenes > 64)
				numGenes = 64;
			if(numGenes < 1)
				numGenes = 1;
			
			HashMap<Integer, Gene> genome = new HashMap<Integer, Gene>();
			
			for(int j = 0; j < numGenes; j++){
				
				int next = random.nextInt(128);
				Gene gene = geneList.get(next);
				
				//if the gene is even -- positive -- and the genome contains neither it nor its negative
				if(next % 2 == 0 && !(genome.containsKey(next) || genome.containsKey(next + 1))){

					gene.setWeight((random.nextDouble() + 0.1)*gene.getWeight());
					genome.put(next, gene);
					
				}
				else if(next % 2 != 0 && !(genome.containsKey(next) || genome.containsKey(next - 1))){

					gene.setWeight((random.nextDouble() + 0.1)*gene.getWeight());
					genome.put(next, gene);
					
				}
				else
					j--;
				
			}
			
			generation.put(i, genome);
			
		}
		
		double[][] distanceMap = distanceMapFirstGen();
		speciateFirstGen(distanceMap);
		fitness(distanceMap);
		
		System.out.println(delta + ":" + speciesNum);
		
	}
	
	private double[][] distanceMapFirstGen(){
		
		double[][] distanceMap = new double[num][num];
		
		for(int i = 0; i < num; i++){
			
			for(int j = 0; j < num; j++){
				
				if(i != j && distanceMap[j][i] == 0){
					
					HashMap<Integer, Gene> genomeI = generation.get(i);
					HashMap<Integer, Gene> genomeJ = generation.get(j);
					
					int N = Math.max(genomeI.size(), genomeJ.size());
					
					Object[] keysI = genomeI.keySet().toArray();
					Object[] keysJ = genomeJ.keySet().toArray();
					
					double Wn = 0;
					double n = 0;
					double D = 0;
					double E = 0;
					double maxI = 0;
					double maxJ = 0;
					
					for(int a = 0; a < genomeI.size(); a++){
						
						for(int b = 0; b < genomeJ.size(); b++){
							
							maxI = (double) Math.max((Integer)keysI[a], maxI);
							maxJ = (double) Math.max((Integer)keysJ[b], maxJ);
							
							if((Integer)keysI[a] == (Integer)keysJ[b]){
								
								n++;
								Wn += Math.abs(genomeI.get((Integer)keysI[a]).getWeight() - genomeJ.get((Integer)keysJ[b]).getWeight());
								
							}
							
						}
						
					}
					
					E = Math.abs(maxI - maxJ);
					D = Math.min(maxI, maxJ) - n;
					
					double W;
					
					if(n == 0)
						W = Integer.MAX_VALUE;
					else
						W = Wn/n;
					
					distanceMap[i][j] = E/N + D/N + W;
					distanceMap[j][i] = distanceMap[i][j];
					
				}
				
			}
			
		}
		
		PriorityQueue<Double> queue = new PriorityQueue<Double>();
		
		for(int i = 0; i < distanceMap.length; i++){
			
			for(int j = 0; j < distanceMap[i].length; j++){
				
				queue.add(distanceMap[i][j]);
				
			}
			
		}
		
		int size = queue.size();
		for(int i = 0; i < size; i++){
		
			if(i == 1500)
				delta = queue.poll();
			else
				queue.poll();
			
		}
		
		return distanceMap;
		
	}

	private void speciateFirstGen(double[][] distanceMap){
		
		species.put(0, 0);
		speciesNum++;
		
		for(int i = 1; i < distanceMap.length; i++){
			
			for(int j = 0; j < i; j++){
				
				if(distanceMap[i][j] <= delta){
					
					species.put(i, species.get(j));
					break;
					
				}
				else if(j == i - 1)
					species.put(i, speciesNum++);
				
			}
			
		}
		
	}
	
	private void fitness(double[][] distanceMap){
		
		int n = 10;
		
		for(int i = 0; i < num; i++){
			
			Brain brain = genoToPheno(generation.get(i));
			
			for(int j = 0; j < n; j++){
				
				int score = 0;
				
				while(game.canWin){
					
					int[][] boardB4 = game.getBoard().clone();
					
					score = game.sumAll();
					
					game.move(brain.brainMove(game.getBoard()));
					
					//If it doesn't make a move that changes anything
					if(boardB4.equals(game.getBoard()))
						break;
					
				}
				if(fitness.containsKey(i))
					score += fitness.get(i);
				fitness.put(i, (double)score);
				game.resetBoard();
				
			}
			
			int neighbors = 0;
			
			for(int j = 0; j < num; j++){
				
				if(distanceMap[i][j] <= delta)
					neighbors++;
				
			}
			
			fitness.put(i, fitness.get(i)/n);
			//fitness.put(i, fitness.get(i)/neighbors);
			System.out.println(fitness.get(i));
			
		}
		
	}
	
	private Brain genoToPheno(HashMap<Integer, Gene> genome){
		
		double[][] graph = new double[nodes][nodes];
		
		Gene[] genes = genome.values().toArray(new Gene[0]);
		
		for(int i = 0; i < genes.length; i++){
			
			graph[genes[i].getIn()][genes[i].getOut()] = genes[i].getWeight();
			
		}
		
		return new Brain(nodes, graph);
		
	}
	
}