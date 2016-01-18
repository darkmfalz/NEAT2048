import java.util.*;

public class Breeder{
	
	private int num;
	private int speciesNum;
	private int nodes;
	private int innovTotal;
	private int generations;
	private double delta = 10;
	private HashMap<Integer, Gene> geneList;
	private HashMap<Integer, HashMap<Integer, Gene>> generation;
	private HashMap<Integer, Double> tiers;
	private int[] species; //The key is the organism number, the value is the species number
	private double[] fitness; //The key is organism number, the value is fitness number
	private Game game;
	
	public Breeder(int num){
		
		this.num = num;
		speciesNum = 0;
		nodes = 20;
		innovTotal = 0;
		generations = 0;
		
		geneList = new HashMap<Integer, Gene>();
		generation = new HashMap<Integer, HashMap<Integer, Gene>>();
		tiers = new HashMap<Integer, Double>();
		species = new int[num];
		fitness = new double[num];
		
		game = new Game();
		
		makeFirstGen();
		
		for(int i = 0; i < 100; i++)
			breedNextGen();
		
	}
	
	private void breedNextGen(){
		
		int organisms = 0;
		
		generations++;
		
		HashMap<Integer, HashMap<Integer, Gene>> prevGeneration = new HashMap<Integer, HashMap<Integer, Gene>>();
		generation.putAll(prevGeneration);
		generation.clear();
		
		for(int i = 0; i < speciesNum; i++){
			
			PriorityQueue<Double> thisSpecies = new PriorityQueue<Double>();
			
			for(int j = 0; j < species.length; j++){
				
				if(species[j] == i)
					thisSpecies.add(fitness[j]);
				
			}
			
			int[] thisSpeciesArr = new int[thisSpecies.size()/2];
			
			for(int j = 0; j < thisSpecies.size()/2; j++){
				
				double currFitness = thisSpecies.poll();
				
				for(int k = 0; k < fitness.length; k++){
					
					if(species[k] == i && fitness[k] == currFitness)
						thisSpeciesArr[j] = k;
					
				}
				
			}
			
			Random random = new Random();
			
			for(int j = 0; j < thisSpeciesArr.length; j++){
				
				int parentID1 = random.nextInt(thisSpeciesArr.length - 1);
				if(parentID1 >= j)
					parentID1++;
				int parentID2 = random.nextInt(thisSpeciesArr.length - 1);
				if(parentID2 >= j)
					parentID2++;
				
				HashMap<Integer, Gene> parent = prevGeneration.get(thisSpeciesArr[j]);
				HashMap<Integer, Gene> parent1 = prevGeneration.get(parentID1);
				HashMap<Integer, Gene> parent2 = prevGeneration.get(parentID2);
				
				HashMap<Integer, Gene> child1 = makeChild(parent, parent1);
				HashMap<Integer, Gene> child2 = makeChild(parent, parent2);
			
				generation.put(organisms++, child1);
				generation.put(organisms++, child2);
				
			}
			
		}
		
		double[][] distanceMap = distanceMap(prevGeneration);
		System.out.println("Generation:" + generations + ", Delta:" + delta + ", Species:" + speciesNum);
		fitness(distanceMap);
		
	}
	
	private void makeFirstGen(){
		
		for(int i = 0; i < 4; i++)
			tiers.put(i, ((double) Integer.MAX_VALUE)/2);
		for(int i = 4; i < 20; i++)
			tiers.put(i, (double) 0);
		
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
		System.out.println("Generation:" + generations + ", Delta:" + delta + ", Species:" + speciesNum);
		fitness(distanceMap);
		
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
		
		species[0] = 0;
		speciesNum++;
		
		for(int i = 1; i < distanceMap.length; i++){
			
			for(int j = 0; j < i; j++){
				
				if(distanceMap[i][j] <= delta){
					
					species[i] = species[j];
					break;
					
				}
				else if(j == i - 1)
					species[i] = speciesNum++;
				
			}
			
		}
		
	}
	
	private double[][] distanceMap(HashMap<Integer, HashMap<Integer, Gene>> prevGeneration){
		
		double[][] distanceMap = new double[num][num];
		
		for(int i = 0; i < num; i++){
			
			for(int j = 0; j < num; j++){
				
				HashMap<Integer, Gene> genomeI = generation.get(i);
				HashMap<Integer, Gene> genomeJ = prevGeneration.get(j);
				
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
	
	private void fitness(double[][] distanceMap){
		
		double maxfit = 0;
		int maxfitspec = 0;
		double minfit = Double.MAX_VALUE;
		int minfitspec = 0;
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
				score += fitness[i];
				fitness[i] = (double)score;
				game.resetBoard();
				
			}
			
			int neighbors = 0;
			
			for(int j = 0; j < num; j++){
				
				if(distanceMap[i][j] <= delta)
					neighbors++;
				
			}
			
			fitness[i] = fitness[i]/n;
			//fitness[i] = fitness[i]/neighbors;
			
			if(fitness[i] > maxfit){
				
				maxfit = fitness[i];
				maxfitspec = i;
				
			}
			
			if(fitness[i] < minfit){
				
				minfit = fitness[i];
				minfitspec = i;
				
			}
			
		}
		
		System.out.println("Organism " + maxfitspec + "\nfitness:" + fitness[maxfitspec] + "\ntopology:");
		double[][] fit = genoToPheno(generation.get(maxfitspec)).getGraph().clone();
		
		System.out.print("{");
		for(int i = 0; i < fit.length; i++){
			
			System.out.print("{");
			for(int j = 0; j < fit[i].length; j++){
				
				if(j != fit[i].length - 1)
					System.out.print(fit[i][j] + ",");
				else
					System.out.print(fit[i][j]);
				
			}
			if(i != fit.length - 1)
				System.out.print("},\n");
			else
				System.out.print("}");
			
		}
		System.out.println("}");
		
		//Prints out the LEAST fit organism
		/*System.out.println("Organism " + minfitspec + "\nfitness:" + fitness[minfitspec] + "\ntopology:");
		fit = genoToPheno(generation.get(minfitspec)).getGraph().clone();
		
		System.out.print("{");
		for(int i = 0; i < fit.length; i++){
			
			System.out.print("{");
			for(int j = 0; j < fit[i].length; j++){
				
				if(j != fit[i].length - 1)
					System.out.print(fit[i][j] + ",");
				else
					System.out.print(fit[i][j]);
				
			}
			if(i != fit.length - 1)
				System.out.print("},\n");
			else
				System.out.print("}");
			
		}
		System.out.println("}");*/
		
	}
	
	private Brain genoToPheno(HashMap<Integer, Gene> genome){
		
		double[][] graph = new double[nodes][nodes];
		
		Gene[] genes = genome.values().toArray(new Gene[0]);
		
		for(int i = 0; i < genes.length; i++){
			
			graph[genes[i].getIn()][genes[i].getOut()] = genes[i].getWeight();
			
		}
		
		return new Brain(nodes, graph);
		
	}
	
	private HashMap<Integer, Gene> makeChild(HashMap<Integer, Gene> parent1, HashMap<Integer, Gene> parent2){
		
		Random random = new Random();
		
		HashMap<Integer, Gene> child = new HashMap<Integer, Gene>();
		
		Integer[] parentKeys1 = parent1.keySet().toArray(new Integer[0]);
		Integer[] parentKeys2 = parent2.keySet().toArray(new Integer[0]);
		
		for(int i = 0; i < parentKeys1.length; i++){
			
			for(int j = 0; j < parentKeys2.length; j++){
				
				if(parentKeys1[i] == parentKeys2[j]){
					
					boolean rand = random.nextBoolean();
					
					if(rand)
						child.put(parentKeys1[i], parent1.get(i));
					else
						child.put(parentKeys2[j], parent2.get(j));
					
				}
				
			}
			
		}
		
		double mutateDrawing = random.nextDouble();
		
		if(mutateDrawing <= 0.25){
		
			int n = (int) Math.round(Math.abs(random.nextGaussian()) + 1);
			
			for(int i = 0; i < n; i++)
				child = mutate(child);
			
		}
		
		return child;
		
	}

	private HashMap<Integer, Gene> mutate(HashMap<Integer, Gene> organism){
		
		Random random = new Random();
		
		boolean mutationType = random.nextBoolean();
		
		//Add Connection
		if(mutationType){
			
			Integer[] keys = organism.keySet().toArray(new Integer[0]);
			
			HashMap<Integer, Double> nodes = new HashMap<Integer, Double>();
			for(int i = 0; i < keys.length; i++){
				
				if(!nodes.containsKey(keys[i]))
					nodes.put(keys[i], tiers.get(keys[i]));
				
			}
			
			keys = nodes.keySet().toArray(new Integer[0]);
			int in = -1;
			int out = -1;
			boolean cont = true;
			while(cont){
				
				in = keys[random.nextInt(keys.length)];
				out = keys[random.nextInt(keys.length)];
				
				if(nodes.get(in) > nodes.get(out)){
					
					cont = false;
					
				}
				
			}
			
			boolean doesntExist = true;
			Gene newPosGene = new Gene(in, out, 1.0, innovTotal);
			Gene[] genes = geneList.values().toArray(new Gene[0]);
			for(int i = 0; i < genes.length; i++){
				
				if(newPosGene.getIn() == genes[i].getIn() && newPosGene.getOut() == genes[i].getOut()){
					
					if(random.nextBoolean())
						organism.put(genes[i].getInnov(), new Gene(in, out, random.nextDouble() + .1, genes[i].getInnov()));
					else
						organism.put(genes[i].getInnov() + 1, new Gene(in, out, -1*random.nextDouble() - .1, genes[i].getInnov() + 1));
					doesntExist = false;
					
				}
					
			}
			
			if(doesntExist){
				
				geneList.put(innovTotal++, newPosGene);
				Gene newNegGene = new Gene(in, out, -1.0, innovTotal);
				geneList.put(innovTotal++, newNegGene);
				
				if(random.nextBoolean())
					organism.put(newPosGene.getInnov(), new Gene(in, out, random.nextDouble() + .1, newPosGene.getInnov()));
				else
					organism.put(newNegGene.getInnov(), new Gene(in, out, -1*random.nextDouble() - .1, newNegGene.getInnov()));
				
			}
			
		}
		//Add Node
		else{
			
			Integer[] keys = organism.keySet().toArray(new Integer[0]);
			int genePos = random.nextInt(keys.length);
			int geneID = keys[genePos];
			Gene gene = organism.get(geneID);
			organism.remove(geneID);
			
			double tierOut = tiers.get(gene.getOut());
			double tierIn = tiers.get(gene.getIn());
			
			double tier = (tierOut + tierIn)/2;
			
			int in = gene.getIn();
			int out = gene.getOut();
			
			Gene newPosGene1 = new Gene(nodes, out, 1.0, innovTotal);
			geneList.put(innovTotal++, newPosGene1);
			Gene newNegGene1 = new Gene(nodes, out, -1.0, innovTotal);
			geneList.put(innovTotal++, newNegGene1);
			Gene newPosGene2 = new Gene(in, nodes, 1.0, innovTotal);
			geneList.put(innovTotal++, newPosGene2);
			Gene newNegGene2 = new Gene(in, nodes, -1.0, innovTotal);
			geneList.put(innovTotal++, newNegGene2);
			
			tiers.put(nodes++, tier);
			
			if(gene.getWeight() > 0.0){
				
				Gene gene1 = new Gene(nodes, out, 1.0, newPosGene1.getInnov());
				gene1.setWeight(gene.getWeight());
				Gene gene2 = new Gene(in, nodes, 1.0, newPosGene2.getInnov());
				gene2.setWeight(gene.getWeight());
				
				organism.put(gene1.getInnov(), gene1);
				organism.put(gene2.getInnov(), gene2);
				
			}
			else{
				
				Gene gene1 = new Gene(nodes, out, -1.0, newNegGene1.getInnov());
				gene1.setWeight(gene.getWeight());
				Gene gene2 = new Gene(in, nodes, -1.0, newNegGene2.getInnov());
				gene2.setWeight(gene.getWeight());
				
				organism.put(gene1.getInnov(), gene1);
				organism.put(gene2.getInnov(), gene2);
				
			}
			
		}
		
		return organism;
		
	}
	
}