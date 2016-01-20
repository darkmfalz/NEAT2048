import java.util.*;

public class Breeder{

	private int numOrgo;
	private int oldNumOrgo;
	private int speciesNum;
	private int nodes;
	private int innovTotal;
	private int generations;
	public double maxFitness = 0;;
	private double delta = 10;
	private HashMap<Integer, Gene> geneList;
	private HashMap<Integer, Organism> generation;
	private HashMap<Integer, Double> tiers;
	private Game game;
	
	public Breeder(int numOrgo){
		
		this.numOrgo = numOrgo;
		oldNumOrgo = numOrgo;
		speciesNum = 0;
		nodes = 20;
		innovTotal = 0;
		generations = 0;
		
		geneList = new HashMap<Integer, Gene>();
		generation = new HashMap<Integer, Organism>();
		tiers = new HashMap<Integer, Double>();
		
		game = new Game();
		
		makeFirstGen();
		
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
		
		for(int i = 0; i < numOrgo; i++){
			
			int numGenes = 32 + (int)(16.0*random.nextGaussian());
			
			if(numGenes > 64)
				numGenes = 64;
			if(numGenes < 1)
				numGenes = 1;
			
			HashMap<Integer, Gene> genome = new HashMap<Integer, Gene>();
			
			for(int j = 0; j < numGenes; j++){
				
				int next = random.nextInt(128);
				Gene gene = new Gene(geneList.get(next).getIn(), geneList.get(next).getOut(), geneList.get(next).getWeight(), geneList.get(next).getInnov());
				
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
			
			generation.put(i, new Organism(i, genome));
			
		}
		
		double[][] distanceMap = distanceMapFirstGen();
		speciateFirstGen(distanceMap);
		System.out.println("Generation:" + generations + ", Delta:" + delta + ", Species:" + speciesNum + ", Growth:" + (numOrgo - oldNumOrgo));
		oldNumOrgo = numOrgo;
		fitness();
		
	}
	
	private double[][] distanceMapFirstGen(){
		
		double[][] distanceMap = new double[numOrgo][numOrgo];
		
		for(int i = 0; i < numOrgo; i++){
			
			for(int j = 0; j < numOrgo; j++){
				
				if(i != j && distanceMap[j][i] == 0){
					
					HashMap<Integer, Gene> genomeI = generation.get(i).cloneGenome();
					HashMap<Integer, Gene> genomeJ = generation.get(j).cloneGenome();
					
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
		
			double limit = 1500.0/10000.0*numOrgo*oldNumOrgo;
			if(i == (int)limit)
				delta = queue.poll();
			else
				queue.poll();
			
		}
		
		return distanceMap;
		
	}

	private void speciateFirstGen(double[][] distanceMap){
		
		generation.get(0).setSpecies(0);
		speciesNum++;
		
		for(int i = 1; i < distanceMap.length; i++){
			
			int[] speciated = new int[i];
			for(int j = 0; j < i; j++)
				speciated[j] = j;
			
			speciated = shuffleArray(speciated);
			
			for(int j = 0; j < i; j++){
				
				int comp = speciated[j];
				
				if(distanceMap[i][comp] <= delta){
					
					generation.get(i).setSpecies(generation.get(comp).getSpecies());
					break;
					
				}
				else if(j == i - 1)
					generation.get(i).setSpecies(speciesNum++);
				
			}
			
		}
		
		int[] speciesN = new int[speciesNum];
		for(int i = 0; i < generation.size(); i++)
			speciesN[generation.get(i).getSpecies()]++;
		Random random = new Random();
		for(int i = 0; i < speciesN.length; i++){
			
			//if there are an odd number of organisms
			//randomly clone
			if(speciesN[i] % 2 == 1){
				
				int randID = random.nextInt(generation.size());
				while(generation.get(randID).getSpecies() != i)
					randID = random.nextInt(generation.size() - 1);
				Organism twin = generation.get(randID).clone();
				twin.setID(numOrgo);
				generation.put(numOrgo, twin);
				
				numOrgo++;
				speciesN[i]++;
				
			}
			
		}
		/*
		//eliminate species w/ < 4 organisms
		for(int i = 0; i < speciesN.length; i++){
			
			if(speciesN[i] < 4){
				
				//remove organisms
				for(int j = 0; j < generation.size(); j++){
					
					if(generation.get(j).getSpecies() == i){
						
						generation.remove(j);
						speciesN[i]--;
						//replace with random twin
						int randID = random.nextInt(generation.size() - 1);
						while(generation.get(randID).getSpecies() == i)
							randID = random.nextInt(generation.size() - 1);
						if(randID >= j)
							randID++;
						Organism twin = generation.get(randID).clone();
						twin.setID(j);
						generation.put(j, twin);
						speciesN[twin.getSpecies()]++;
						
					}
					
				}
				
			}
			
		}*/
		//Or, instead, add the necessary number to continue the existence of species
		for(int i = 0; i < speciesN.length; i++){
			
			//if there are an odd number of organisms
			//randomly clone
			while(speciesN[i] < 4){
				
				int randID = random.nextInt(generation.size());
				while(generation.get(randID).getSpecies() != i)
					randID = random.nextInt(generation.size() - 1);
				Organism twin = generation.get(randID).clone();
				twin.setID(numOrgo);
				generation.put(numOrgo, twin);
				
				numOrgo++;
				speciesN[i]++;
				
			}
			
		}
		
	}
	
	public void breedNextGen(){
		
		int organisms = 0;
		
		HashMap<Integer, Organism> prevGeneration = new HashMap<Integer, Organism>(generation);
		oldNumOrgo = prevGeneration.size();
		generation.clear();
		
		for(int i = 0; i < speciesNum; i++){

			Comparator<Organism> fitnessComparator = new Comparator<Organism>(){
				@Override
				public int compare(Organism o1, Organism o2){
					return -1*Double.compare(o1.getFitness(), o2.getFitness());
					}
				};
				
			PriorityQueue<Organism> fitnessQueue = new PriorityQueue<Organism>(fitnessComparator);
			
			for(int j = 0; j < prevGeneration.size(); j++)
				if(prevGeneration.get(j).getSpecies() == i)
					fitnessQueue.add(prevGeneration.get(j));

			Organism[] fit = new Organism[fitnessQueue.size()/2];
			for(int j = 0; j < fit.length; j++)
				fit[j] = fitnessQueue.poll();
			
			Random random = new Random();
			for(int j = 0; j < fit.length; j++){
				
				Organism parent = fit[j];
				
				int parentID1 = random.nextInt(fit.length - 1);
				if(parentID1 >= j)
					parentID1++;
				Organism parent1 = fit[parentID1];
				
				int parentID2 = random.nextInt(fit.length - 1);
				if(parentID2 >= j)
					parentID2++;
				Organism parent2 = fit[parentID2];
				
				Organism child1 = makeChild(parent, parent1);
				child1.setID(organisms);
				generation.put(organisms, child1);
				organisms++;
				
				Organism child2 = makeChild(parent, parent2);
				child2.setID(organisms);
				generation.put(organisms, child2);
				organisms++;
		
			}
			
		}
		
		generations++;
		double[][] distanceMap = distanceMap(prevGeneration);
		speciate(distanceMap);
		System.out.println("Generation:" + generations + ", Delta:" + delta + ", Species:" + speciesNum + ", Growth:" + (numOrgo - oldNumOrgo));
		fitness();
		
	}

	private double[][] distanceMap(HashMap<Integer, Organism> prevGeneration){
		
		double[][] distanceMap = new double[numOrgo][oldNumOrgo];
		
		for(int i = 0; i < numOrgo; i++){
			
			for(int j = 0; j < oldNumOrgo; j++){
			
				HashMap<Integer, Gene> genomeI = generation.get(i).cloneGenome();
				HashMap<Integer, Gene> genomeJ = prevGeneration.get(j).cloneGenome();
				
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
		
			double limit = 1500.0/10000.0*numOrgo*oldNumOrgo;
			if(i == (int)limit)
				delta = queue.poll();
			else
				queue.poll();
			
		}
		
		return distanceMap;
		
	}

	private void speciate(double[][] distanceMap){
		
		for(int i = 0; i < distanceMap.length; i++){
			
			int[] speciated = new int[distanceMap[i].length];
			for(int j = 0; j < distanceMap[i].length; j++)
				speciated[j] = j;
			
			speciated = shuffleArray(speciated);
			
			for(int j = 0; j < distanceMap[i].length; j++){
				
				int comp = speciated[j];
				
				if(distanceMap[i][comp] <= delta){
					
					generation.get(i).setSpecies(generation.get(comp).getSpecies());
					break;
					
				}
				else if(j == i - 1)
					generation.get(i).setSpecies(speciesNum++);
				
			}
			
		}
		
		int[] speciesN = new int[speciesNum];
		for(int i = 0; i < generation.size(); i++)
			speciesN[generation.get(i).getSpecies()]++;
		Random random = new Random();
		for(int i = 0; i < speciesN.length; i++){
			
			//if there are an odd number of organisms
			//randomly clone
			if(speciesN[i] % 2 == 1){
				
				int randID = random.nextInt(generation.size());
				while(generation.get(randID).getSpecies() != i)
					randID = random.nextInt(generation.size() - 1);
				Organism twin = generation.get(randID).clone();
				twin.setID(numOrgo);
				generation.put(numOrgo, twin);
				
				numOrgo++;
				speciesN[i]++;
				
			}
			
		}
		/*
		//eliminate species w/ < 4 organisms
		for(int i = 0; i < speciesN.length; i++){
			
			if(speciesN[i] < 4){
				
				//remove organisms
				for(int j = 0; j < generation.size(); j++){
					
					if(generation.get(j).getSpecies() == i){
						
						generation.remove(j);
						speciesN[i]--;
						//replace with random twin
						int randID = random.nextInt(generation.size() - 1);
						while(generation.get(randID).getSpecies() == i)
							randID = random.nextInt(generation.size() - 1);
						if(randID >= j)
							randID++;
						Organism twin = generation.get(randID).clone();
						twin.setID(j);
						generation.put(j, twin);
						speciesN[twin.getSpecies()]++;
						
					}
					
				}
				
			}
			
		}*/
		//Or, instead, add the necessary number to continue the existence of species
		for(int i = 0; i < speciesN.length; i++){
			
			//if there are an odd number of organisms
			//randomly clone
			while(speciesN[i] < 4 && speciesN[i] > 0){
				
				int[] speciesArr = new int[speciesN[i]];
				int counter = 0;
				for(int j = 0; j < generation.size(); j++)
					if(generation.get(j).getSpecies() == i)
						speciesArr[counter++] = j;	
				speciesArr = shuffleArray(speciesArr);
				int randID = speciesArr[0];
				
				Organism twin = generation.get(randID).clone();
				twin.setID(numOrgo);
				generation.put(numOrgo, twin);
				
				numOrgo++;
				speciesN[i]++;
				
			}
			
		}
		
	}
	
	private void fitness(){
		
		double maxfit = 0;
		int maxfitspec = 0;
		double minfit = Double.MAX_VALUE;
		int minfitspec = 0;
		int n = 100;
		
		for(int i = 0; i < numOrgo; i++){
			
			Brain brain = genoToPheno(generation.get(i).cloneGenome());
			
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
				score += generation.get(i).getFitness();
				generation.get(i).setFitness((double)score);
				game.resetBoard();
				
			}
			
			/**int neighbors = 0;
			
			for(int j = 0; j < oldNumOrgo; j++){
				
				if(distanceMap[i][j] <= delta)
					neighbors++;
				
			}**/
			
			generation.get(i).setFitness(generation.get(i).getFitness()/n);
			//fitness[i] = fitness[i]/neighbors;
			
			if(generation.get(i).getFitness() > maxfit){
				
				maxfit = generation.get(i).getFitness();
				maxfitspec = i;
				
			}
			
			if(generation.get(i).getFitness() < minfit){
				
				minfit = generation.get(i).getFitness();
				minfitspec = i;
				
			}
			
		}
		
		maxFitness = generation.get(maxfitspec).getFitness();
		
		System.out.println("Organism " + maxfitspec + ", Species " + generation.get(maxfitspec).getSpecies() + "\nfitness:" + generation.get(maxfitspec).getFitness() + "\ntopology:");
		double[][] fit = genoToPheno(generation.get(maxfitspec).cloneGenome()).getGraph().clone();
		
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
	
	private Organism makeChild(Organism parent1, Organism parent2){
		
		return parent1.clone();
		
	}
	
	private Brain genoToPheno(HashMap<Integer, Gene> genome){
		
		double[][] graph = new double[nodes][nodes];
		
		Gene[] genes = genome.values().toArray(new Gene[0]);
		
		for(int i = 0; i < genes.length; i++){
			
			graph[genes[i].getIn()][genes[i].getOut()] = genes[i].getWeight();
			
		}
		
		return new Brain(nodes, graph);
		
	}
	
	public int[] shuffleArray(int[] ar){
		
		Random rnd = new Random();
		
		for (int i = ar.length - 1; i > 0; i--){
			
			int index = rnd.nextInt(i + 1);
			// Simple swap
			int a = ar[index];
			ar[index] = ar[i];
			ar[i] = a;
			
		}
		
		return ar;
		
	}
	
}