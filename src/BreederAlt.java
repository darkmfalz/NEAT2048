import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.PriorityQueue;
import java.util.Random;

public class BreederAlt{

	public int numOrgo;
	private int maxOrgo = 2000;
	private int oldNumOrgo;
	private int speciesNum;
	private int nodes;
	private int innovTotal;
	public int generations;
	public double maxFitness = 0;
	public double avgFitness = 0;
	public int maxNodes = 0;
	private double delta = 10;
	private HashMap<Integer, Gene> geneList;
	private HashMap<Integer, Organism> generation;
	private HashMap<Integer, Double> tiers;
	private Game game;
	private boolean oneChildPolicy;
	
	private double proprSpecies = 0.25;
	
	public BreederAlt(int numOrgo){
		
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
		
		oneChildPolicy = false;
		
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
		System.out.println("Generation:" + generations + ", Delta:" + delta + ", Organisms: " + numOrgo + ", Species:" + speciesNum + ", Growth:" + (numOrgo - oldNumOrgo));oldNumOrgo = numOrgo;
		distanceMap = distanceMapFirstGen();
		fitness(distanceMap);
		
		if(generation.size() > maxOrgo)
			oneChildPolicy = true;
		
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
		
		//Need better delta selection
		PriorityQueue<Double> queue = new PriorityQueue<Double>();
		PriorityQueue<Double> queue2 = new PriorityQueue<Double>();
		
		for(int i = 0; i < distanceMap.length; i++){
			
			for(int j = 0; j < distanceMap[i].length; j++){
				
				queue.add(distanceMap[i][j]);
				
			}
			int size = queue.size();
			for(int j = 0; j < size; j++){
				
				if(j < size*proprSpecies+1 && j > size*proprSpecies-1)
					queue2.add(queue.poll());
				else
					queue.poll();
				
			}
			
		}
		
		int size = queue2.size();
		for(int i = 0; i < size; i++){
	
			if(i == (int)size/2)
				delta = queue2.poll();
			else
				queue2.poll();
			
		}
		
		delta = Math.min(delta, Integer.MAX_VALUE/2);
		
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
		
		Comparator<Organism> fitnessComparator = new Comparator<Organism>(){
			@Override
			public int compare(Organism o1, Organism o2){
				return -1*Double.compare(o1.getFitness(), o2.getFitness());
				}
			};
			
		PriorityQueue<Organism> fitnessQueue = new PriorityQueue<Organism>(fitnessComparator);
		
		for(int i = 0; i < prevGeneration.size(); i++)
			fitnessQueue.add(prevGeneration.get(i));
		
		Organism[] orgo = new Organism[prevGeneration.size()/2 + 1];
		for(int i = 0; i < orgo.length; i++)
			orgo[i] = fitnessQueue.poll();
		
		Random random = new Random();
		for(int i = 0; i < orgo.length; i++){
			
			ArrayList<Organism> species = new ArrayList<Organism>();
			for(int j = 0; j < orgo.length; j++)
				if(orgo[i].getSpecies() == orgo[j].getSpecies())
					species.add(orgo[j]);
					
			if(species.size() > 1){
				
				Organism parent = orgo[i];
				
				int parentID1 = random.nextInt(species.size());
				while(species.get(parentID1).getID() == orgo[i].getSpecies())
					parentID1 = random.nextInt(species.size());
				Organism parent1 = species.get(parentID1);
				
				int parentID2 = random.nextInt(species.size());
				while(species.get(parentID2).getID() == orgo[i].getSpecies())
					parentID2 = random.nextInt(species.size());
				Organism parent2 = species.get(parentID2);
				
				Organism child1 = makeChild(parent, parent1);
				child1.setID(organisms);
				generation.put(organisms, child1);
				organisms++;
				
				if(!oneChildPolicy){
					
					Organism child2 = makeChild(parent, parent2);
					child2.setID(organisms);
					generation.put(organisms, child2);
					organisms++;
					
				}
				
			}
			else{
				
				Organism parent = orgo[i];
				
				Organism child1 = makeChild(parent, parent.clone());
				child1.setID(organisms);
				generation.put(organisms, child1);
				organisms++;
				
			}
				
		}
		
		numOrgo = organisms;
		generations++;
		double[][] distanceMap = distanceMap(prevGeneration);
		speciate(distanceMap, prevGeneration);
		System.out.println("Generation:" + generations + ", Delta:" + delta + ", Organisms: " + numOrgo + ", Species:" + speciesNum + ", Growth:" + (numOrgo - oldNumOrgo));
		fitness(distanceMapFirstGen());
		
		if(generation.size() > maxOrgo)
			oneChildPolicy = true;
		else
			oneChildPolicy = false;
		
	}

	private double[][] distanceMap(HashMap<Integer, Organism> prevGeneration){
		
		double[][] distanceMap = new double[numOrgo][oldNumOrgo];
		
		for(int i = 0; i < generation.size(); i++){
			
			for(int j = 0; j < prevGeneration.size(); j++){
			
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
		
		//Need better delta selection
		PriorityQueue<Double> queue = new PriorityQueue<Double>();
		PriorityQueue<Double> queue2 = new PriorityQueue<Double>();
		
		for(int i = 0; i < distanceMap.length; i++){
			
			for(int j = 0; j < distanceMap[i].length; j++){
				
				queue.add(distanceMap[i][j]);
				
			}
			int size = queue.size();
			for(int j = 0; j < size; j++){
				
				if(j < size*proprSpecies+1 && j > size*proprSpecies-1)
					queue2.add(queue.poll());
				else
					queue.poll();
				
			}
			
		}
		
		int size = queue2.size();
		for(int i = 0; i < size; i++){
	
			if(i == (int)size/2)
				delta = queue2.poll();
			else
				queue2.poll();
			
		}
		
		delta = Math.min(delta, Integer.MAX_VALUE/2);
		
		return distanceMap;
		
	}

	private void speciate(double[][] distanceMap, HashMap<Integer, Organism> prevGeneration){
		
		for(int i = 0; i < distanceMap.length; i++){
			
			int[] speciated = new int[distanceMap[i].length];
			for(int j = 0; j < distanceMap[i].length; j++)
				speciated[j] = j;
			
			speciated = shuffleArray(speciated);
			
			for(int j = 0; j < distanceMap[i].length; j++){
				
				int comp = speciated[j];
				
				if(distanceMap[i][comp] <= delta){
					
					generation.get(i).setSpecies(prevGeneration.get(comp).getSpecies());
					break;
					
				}
				else if(j == distanceMap[i].length - 1)
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
	
	private void fitness(double[][] distanceMap){
		
		avgFitness = 0;
		
		double maxfit = 0;
		int maxfitspec = 0;
		double minfit = Double.MAX_VALUE;
		int minfitspec = 0;
		int n = 30;
		
		for(int i = 0; i < generation.size(); i++){
			
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
			
			generation.get(i).setFitness(generation.get(i).getFitness()/n);
			avgFitness += generation.get(i).getFitness()/generation.size();
			if(generation.get(i).getFitness() > maxfit){
				
				maxfit = generation.get(i).getFitness();
				maxfitspec = i;
				
			}
			
			if(generation.get(i).getFitness() < minfit){
				
				minfit = generation.get(i).getFitness();
				minfitspec = i;
				
			}

			int neighbors = 0;
			
			for(int j = 0; j < numOrgo; j++){
				
				if(distanceMap[i][j] <= delta)
					neighbors++;
				
			}
			generation.get(i).setFitness(generation.get(i).getFitness()/neighbors);
			
		}
		
		maxFitness = maxfit;
		
		System.out.println("Max Number of Nodes:" + maxNodes);
		System.out.println("Organism " + maxfitspec + ", Species " + generation.get(maxfitspec).getSpecies() + "\nfitness:" + maxFitness + "\ntopology:");
		double[][] fit = genoToPheno(generation.get(maxfitspec).cloneGenome()).getGraph().clone();
		
		if(maxFitness > 2000){
		
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
			
		}
		
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
		
		HashMap<Integer, Gene> genome1 = parent1.cloneGenome();
		Gene[] genome1values = genome1.values().toArray(new Gene[0]);
		HashMap<Integer, Gene> genome2 = parent2.cloneGenome();
		Gene[] genome2values = genome2.values().toArray(new Gene[0]);
		//int N = Math.max(genome1.size(), genome2.size());
		
		HashMap<Integer, Gene> genome = new HashMap<Integer, Gene>();
		Random random = new Random();
		
		for(int i = 0; i < genome1.size(); i++){
			
			for(int j = 0; j < genome2.size(); j++){
				
				if(genome1values[i].getInnov() == genome2values[j].getInnov()){
					
					boolean choice = random.nextBoolean();
					
					if(choice){
						double enable = random.nextDouble();
						if(enable > 0.25)
							genome1values[i].enable();
						genome.put(genome1values[i].getInnov(), genome1values[i]);
					}
					else{
						double enable = random.nextDouble();
						if(enable > 0.25)
							genome2values[j].enable();
						genome.put(genome2values[j].getInnov(), genome2values[j]);
					}
					break;
					
				}
				else if(parent1.getFitness() >= parent2.getFitness() && j == genome2.size() - 1){
					double enable = random.nextDouble();
					if(enable > 0.25)
						genome1values[i].enable();
					genome.put(genome1values[i].getInnov(), genome1values[i]);
				}
				
			}
			
		}
		
		if(parent1.getFitness() < parent2.getFitness()){
			
			for(int i = 0; i < genome2.size(); i++){
				
				for(int j = 0; j < genome1.size(); j++){
					
					if(genome1values[j] == genome2values[i])
						break;
					else if(j == genome1.size() - 1){
						double enable = random.nextDouble();
						if(enable > 0.25)
							genome2values[i].enable();
						genome.put(genome2values[i].getInnov(), genome2values[i]);
					}
					
				}
				
			}
			
		}
		
		double mutate = random.nextDouble();
		
		if(mutate <= 0.3)
			genome = mutate(genome);
		
		return new Organism(genome);
		
	}
	
	private HashMap<Integer, Gene> mutate(HashMap<Integer, Gene> genome){
		
		Random random = new Random();
		int mutation = random.nextInt(2);
		//mutation = true;
		
		if(mutation == 0){
			
			//Add Node
			
			Gene[] genomeKeys = genome.values().toArray(new Gene[0]);
			
			Gene gene = genomeKeys[random.nextInt(genomeKeys.length)];
			genome.remove(gene.getInnov());
			gene.disable();
			genome.put(gene.getInnov(), gene);
			
			tiers.put(nodes, (tiers.get(gene.getIn()) + tiers.get(gene.getOut()))/2);
			Gene gene1pos = new Gene(nodes, gene.getOut(), 1.0, innovTotal++);
			Gene gene1neg = new Gene(nodes, gene.getOut(), -1.0, innovTotal++);
			Gene gene2pos = new Gene(gene.getIn(), nodes, 1.0, innovTotal++);
			Gene gene2neg = new Gene(gene.getIn(), nodes, -1.0, innovTotal++);
			
			geneList.put(gene1pos.getInnov(), gene1pos);
			geneList.put(gene1neg.getInnov(), gene1neg);
			geneList.put(gene2pos.getInnov(), gene2pos);
			geneList.put(gene2neg.getInnov(), gene2neg);
			
			nodes++;
			
			if(gene.getWeight() >  0){
				
				Gene gene1 = gene1pos.clone();
				gene1.setWeight(gene.getWeight());
				genome.put(gene1.getInnov(), gene1);
				Gene gene2 = gene2pos.clone();
				gene2.setWeight(gene.getWeight());
				genome.put(gene2.getInnov(), gene2);
				
			}
			else{
				
				Gene gene1 = gene1neg.clone();
				gene1.setWeight(gene.getWeight());
				genome.put(gene1.getInnov(), gene1);
				Gene gene2 = gene2neg.clone();
				gene2.setWeight(gene.getWeight());
				genome.put(gene2.getInnov(), gene2);
				
			}
			
		}
		else{
			
			//Add Edge
			
			//Find out which nodes exist
			ArrayList<Integer> nodes = new ArrayList<Integer>();
			Gene[] genomeKeys = genome.values().toArray(new Gene[0]);
			
			for(int i = 0; i < genomeKeys.length; i++){
				
				if(!nodes.contains(genomeKeys[i].getIn()))
					nodes.add(genomeKeys[i].getIn());
				if(!nodes.contains(genomeKeys[i].getOut()))
					nodes.add(genomeKeys[i].getOut());
				
			}
			
			Integer[] nodesArr = nodes.toArray(new Integer[0]);
			
			int a = random.nextInt(nodes.size());
			int b = random.nextInt(nodes.size() - 1);
			if(b >= a)
				b++;
			int in = nodesArr[a];
			int out = nodesArr[b];
			
			while(tiers.get(in) <= tiers.get(out)){
				
				a = random.nextInt(nodes.size());
				b = random.nextInt(nodes.size() - 1);
				if(b >= a)
					b++;
				in = nodesArr[a];
				out = nodesArr[b];
				
			}
			
			Gene[] geneArr = geneList.values().toArray(new Gene[0]);
			Gene genePos = new Gene(in, out, 1.0, innovTotal);
			
			for(int i = 0; i < geneArr.length; i++){
				
				if(geneArr[i].getIn() == in && geneArr[i].getOut() == out && geneArr[i].getWeight() > 0){
					
					genePos = geneArr[i];
					break;
				
				}
				
			}
			
			Gene geneNeg = new Gene(in, out, -1.0, genePos.getInnov() + 1);
			
			if(genePos.getInnov() == innovTotal){
				
				geneList.put(genePos.getInnov(), genePos);
				geneList.put(geneNeg.getInnov(), geneNeg);
				
				innovTotal++;
				innovTotal++;
			
			}
			
			boolean pos = random.nextBoolean();
			
			if(!genome.containsKey(genePos.getInnov()) && !genome.containsKey(geneNeg.getInnov())){
			
				if(pos){
					
					double weight = (random.nextDouble() + .1)*genePos.getInnov();
					Gene gene = new Gene(in, out, weight, genePos.getInnov());
					genome.put(gene.getInnov(), gene);
					
				}
				else{
					
					double weight = (random.nextDouble() + .1)*geneNeg.getInnov();
					Gene gene = new Gene(in, out, weight, geneNeg.getInnov());
					genome.put(gene.getInnov(), gene);
					
				}
				
			}
			else{
				
				if(genome.containsKey(genePos.getInnov()))
					genome.get(genePos.getInnov()).enable();
				else
					genome.get(geneNeg.getInnov()).enable();
				
			}
			
		}
		
		return genome;
		
	}
	
 	private Brain genoToPheno(HashMap<Integer, Gene> genome){
		
 		Gene[] genes = genome.values().toArray(new Gene[0]);
		
		PriorityQueue<Integer> nodesL = new PriorityQueue<Integer>();
		Gene[] genomeKeys = genome.values().toArray(new Gene[0]);
		
		for(int i = 0; i < 20; i++)
			nodesL.add(i);
		for(int i = 20; i < genomeKeys.length; i++){
			
			if(!nodesL.contains(genomeKeys[i].getIn()))
				nodesL.add(genomeKeys[i].getIn());
			if(!nodesL.contains(genomeKeys[i].getOut()))
				nodesL.add(genomeKeys[i].getOut());
			
		}
		
		maxNodes = Math.max(maxNodes, nodesL.size());
		
		double[][] graph = new double[nodesL.size()][nodesL.size()];
		int[] vertices = new int[nodesL.size()];
	
		for(int i = 0; i < vertices.length; i++)
			vertices[i] = nodesL.poll();
		
		for(int i = 0; i < genes.length; i++){
			
			if(genes[i].getEnabled()){
				
				for(int j = 0; j < vertices.length; j++){
					
					if(vertices[j] == genes[i].getIn()){
						
						for(int k = 0; k < vertices.length; k++){
							
							if(vertices[k] == genes[i].getOut()){
								
								graph[j][k] = genes[i].getWeight();
								break;
								
							}
							
						}
						break;
						
					}
					
				}
				
			}
			
		}

		return new Brain(graph.length, graph);
		
	}
	
 	public void printFitness(){
 		
 		for(int i = 0; i < generation.size(); i++)
 			System.out.println(generation.get(i).getFitness());
 		
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
	
	public double[][] removeRow(double[][] array, int row){
		
		double[][] newArray = new double[array.length - 1][];
		
		for(int i = 0; i < array.length; i++){
			
			if(i != row && i < row)
				newArray[i] = array[i];
			else if(i != row)
				newArray[i - 1] = array[i];
			
		}
		
		return newArray;
		
	}

	public double[][] removeColumn(double[][] array, int column){
		
		double[][] newArray = new double[array.length][array[0].length - 1];
		
		for(int i = 0; i < array.length; i++){
			
			for(int j = 0; j < array[i].length; j++){
				
				if(j != column && j < column)
					newArray[i][j] = array[i][j];
				else if(j != column)
					newArray[i][j - 1] = array[i][j];
				
			}
			
		}
		
		return newArray;
		
	}
	
}