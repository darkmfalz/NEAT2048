import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.PriorityQueue;
import java.util.Random;

public class Breeder {
	
	private int innovTotal;
	public int generations;
	public double maxFitness = 0;
	private double delta = 10;
	private HashMap<Integer, Gene> geneList;
	private HashMap<Integer, Organism> generation;
	private HashMap<Integer, Integer> speciesSize;
	private HashMap<Integer, double[]> extinctionMap;
	private int activeSpecies;
	private int nextSpecies;
	private Game game;

	private double proprSpecies = 0.5;
	
	public Breeder(int numOrgo){
		
		innovTotal = 0;
		generations = 0;
		
		geneList = new HashMap<Integer, Gene>();
		generation = new HashMap<Integer, Organism>();
		
		speciesSize = new HashMap<Integer, Integer>();
		activeSpecies = 0;
		nextSpecies = 0;
		
		extinctionMap = new HashMap<Integer, double[]>();
		
		game = new Game();
		
		makeFirstGen(numOrgo);
		
	}
	
	private void makeFirstGen(int numOrgo){
		
		for(int i = 0; i < 4; i++)
			geneList.put(innovTotal, new Gene(innovTotal++, false, -1, -1, 0, i, Double.MAX_VALUE));
		for(int i = 4; i < 20; i++)
			geneList.put(innovTotal, new Gene(innovTotal++, false, -1, -1, 0, i, 0));
		
		//Create all possible starting genes
		for(int i = 0; i < 4; i++){
			
			for(int j = 0; j < 16; j++){
				
				Gene gene = new Gene(innovTotal, true, i, j + 4, 1.0, -1, -1);
				geneList.put(innovTotal, gene);
				innovTotal++;
				
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
				
				int next = 20 + random.nextInt(64);
				while(!geneList.get(next).getIsEdge())
					next = 20 + random.nextInt(64);
				Gene gene = geneList.get(next).clone();
				
				//if the gene is even -- positive -- and the genome contains neither it nor its negative
				if(!genome.containsKey(gene.getInnov())){

					boolean pos = random.nextBoolean();
					
					if(pos)
						gene.setWeight((random.nextDouble() + 0.1));
					else
						gene.setWeight(-1*(random.nextDouble() + 0.1));
					genome.put(next, gene);
					
					Gene neuronIn = geneList.get(gene.getIn()).clone();
					genome.put(neuronIn.getInnov(), neuronIn);
					Gene neuronOut = geneList.get(gene.getOut()).clone();
					genome.put(neuronOut.getInnov(), neuronOut);
					
				}
				else
					j--;
				
			}
			
			generation.put(i, new Organism(i, genome));
			
		}
		
		double[][] distanceMap = distanceMapFirstGen();
		speciateFirstGen(distanceMap);
		System.out.println("Generation:" + generations + ", Delta:" + delta + ", Organisms: " + generation.size() + ", Active Species: " + activeSpecies);
		fitness();
		
	}
	
	private double[][] distanceMapFirstGen(){
		
		double[][] distanceMap = new double[generation.size()][generation.size()];
		
		for(int i = 0; i < generation.size(); i++){
			
			for(int j = 0; j < generation.size(); j++){
				
				if(i != j && distanceMap[j][i] == 0){
					
					HashMap<Integer, Gene> genomeI = generation.get(i).cloneGenome();
					Integer[] keysI = genomeI.keySet().toArray(new Integer[0]);
					HashMap<Integer, Gene> genomeIEdges = new HashMap<Integer, Gene>(genomeI);
					for(int a = 0; a < keysI.length; a++)
						if(!genomeIEdges.get(keysI[a]).getIsEdge())
							genomeIEdges.remove(keysI[a]);
					keysI = genomeIEdges.keySet().toArray(new Integer[0]);
					
					HashMap<Integer, Gene> genomeJ = generation.get(j).cloneGenome();
					Integer[] keysJ = genomeJ.keySet().toArray(new Integer[0]);
					HashMap<Integer, Gene> genomeJEdges = new HashMap<Integer, Gene>(genomeJ);
					for(int b = 0; b < keysJ.length; b++)
						if(!genomeJEdges.get(keysJ[b]).getIsEdge())
							genomeJEdges.remove(keysJ[b]);
					keysJ = genomeJEdges.keySet().toArray(new Integer[0]);
					
					double N = (double) Math.max(keysI.length, keysJ.length);
					
					double Wn = 0;
					double n = 0;
					double D = 0;
					double E = 0;
					double maxI = 0;
					double maxJ = 0;
					
					for(int a = 0; a < genomeIEdges.size(); a++)
						if(genomeIEdges.get(keysI[a]).getIsEdge())
							maxI = (double) Math.max(keysI[a], maxI);
					for(int b = 0; b < genomeJEdges.size(); b++)
						if(genomeJEdges.get(keysJ[b]).getIsEdge())
							maxJ = (double) Math.max(keysJ[b], maxJ);
					
					for(int a = 0; a < keysI.length; a++){
						
						for(int b = 0; b < keysJ.length; b++){
							
							if(keysI[a] == keysJ[b]){
								
								n++;
								Wn += Math.abs(genomeI.get(keysI[a]).getWeight() - genomeJ.get(keysJ[b]).getWeight());
								
							}
							else if(genomeI.get(keysI[a]).getIn() == genomeJ.get(keysJ[b]).getIn() && genomeI.get(keysI[a]).getOut() == genomeJ.get(keysJ[b]).getOut()){
								
								n++;
								Wn += Math.abs(genomeI.get(keysI[a]).getWeight() - genomeJ.get(keysJ[b]).getWeight());
								
							}
							else if(keysI[a] > maxJ && b == 0)
								E++;
							else if(keysJ[b] > maxI && a == 0)
								E++;
							
						}
						
					}
					
					D = Math.min(keysI.length, keysJ.length) - n;
					
					double W;
					
					if(n == 0)
						W = Double.MAX_VALUE;
					else
						W = Wn/n;
					
					distanceMap[i][j] = E/N + D/N + 0.2*W;
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
	
			if(i == (int)(size*0.5))
				delta = queue2.poll();
			else
				queue2.poll();
			
		}
		
		delta = Math.min(delta, 0.8);
		
		return distanceMap;
		
	}

	private void speciateFirstGen(double[][] distanceMap){
		
		activeSpecies = 0;
		speciesSize.clear();
		
		generation.get(0).setSpecies(0);
		speciesSize.put(0, 1);
		activeSpecies++;
		
		for(int i = 1; i < generation.size(); i++){

			HashMap<Integer, Double> avgSpeciesDist = new HashMap<Integer, Double>();
			for(int j = 0; j < i; j++)
				if(avgSpeciesDist.containsKey(generation.get(j).getSpecies()))
					avgSpeciesDist.put(generation.get(j).getSpecies(), avgSpeciesDist.get(generation.get(j).getSpecies()) + distanceMap[i][j]);
				else
					avgSpeciesDist.put(generation.get(j).getSpecies(), distanceMap[i][j]);
			Integer[] keys = avgSpeciesDist.keySet().toArray(new Integer[0]);
			for(int j = 0; j < keys.length; j++)
				if(speciesSize.get(keys[j]) > 0)
					avgSpeciesDist.put(keys[j], avgSpeciesDist.get(keys[j])/(double)speciesSize.get(keys[j]));
			keys = shuffleArray(keys);
			for(int j = 0; j < keys.length; j++){
				
				if(avgSpeciesDist.get(keys[j]) < delta){
					
					generation.get(i).setSpecies(keys[j]);
					speciesSize.put(keys[j], speciesSize.get(keys[j]) + 1);
					break;
					
				}
				else if(j == keys.length - 1){
					
					generation.get(i).setSpecies(speciesSize.size());
					speciesSize.put(speciesSize.size(), 1);
					activeSpecies++;
					
				}
				
			}
			
		}
		
		/*Random random = new Random();
		for(int i = 0; i < speciesSize.size(); i++){

			//if there are an odd number of organisms
			//randomly breed
			if(speciesSize.get(i) % 2 == 1 && speciesSize.get(i) > 0){
				
				int randID1 = random.nextInt(generation.size());
				while(generation.get(randID1).getSpecies() != i)
					randID1 = random.nextInt(generation.size());
				int randID2;
				if(speciesSize.get(i) > 1){
					randID2 = random.nextInt(generation.size());
					while(generation.get(randID2).getSpecies() != i)
						randID2 = random.nextInt(generation.size());
				}
				else
					randID2 = randID1;
				Organism parent1 = generation.get(randID1);
				Organism parent2 = generation.get(randID2);
				Organism child = makeChild(parent1, parent2, false);
				child.setID(generation.size());
				child.setSpecies(i);
				generation.put(generation.size(), child);
				
				speciesSize.put(i, speciesSize.get(i) + 1);
				
			}
			
		}

		for(int i = 0; i < speciesSize.size(); i++){
			
			while(speciesSize.get(i) < 4){
				
				if(speciesSize.get(i) == 0)
					break;
				
				int randID1 = random.nextInt(generation.size());
				while(generation.get(randID1).getSpecies() != i)
					randID1 = random.nextInt(generation.size());
				int randID2;
				if(speciesSize.get(i) > 1){
					randID2 = random.nextInt(generation.size());
					while(generation.get(randID2).getSpecies() != i)
						randID2 = random.nextInt(generation.size());
				}
				else
					randID2 = randID1;
				Organism parent1 = generation.get(randID1);
				Organism parent2 = generation.get(randID2);
				Organism child = makeChild(parent1, parent2, false);
				child.setID(generation.size());
				child.setSpecies(i);
				generation.put(generation.size(), child);
				
				speciesSize.put(i, speciesSize.get(i) + 1);
				
			}
			
		}*/
		
	}
	
	public void breedNextGen(){
		
		HashMap<Integer, Organism> prevGeneration = new HashMap<Integer, Organism>(generation);
		HashMap<Integer, Organism> prevGenerationTemp = new HashMap<Integer, Organism>(generation);
		generation.clear();
		int numOrgo = 0;
		
		Comparator<Organism> fitnessComparator = new Comparator<Organism>(){
			@Override
			public int compare(Organism o1, Organism o2){
				return -1*Double.compare(o1.getFitness(), o2.getFitness());
				}
			};
			
		HashMap<Integer, double[]> fitnessMetrics = new HashMap<Integer, double[]>();
		HashMap<Integer, PriorityQueue<Organism>> c = new HashMap<Integer, PriorityQueue<Organism>>(); 
		
		for(int i = 0; i < prevGeneration.size(); i++)
			if(c.containsKey(prevGeneration.get(i).getSpecies()))
				c.get(prevGeneration.get(i).getSpecies()).add(prevGeneration.get(i));
			else{
				
				PriorityQueue<Organism> queue = new PriorityQueue<Organism>(fitnessComparator);
				queue.add(prevGeneration.get(i));
				c.put(prevGeneration.get(i).getSpecies(), queue);
				
			}
		
		Integer[] keysC = c.keySet().toArray(new Integer[0]);
		for(int i = 0; i < keysC.length; i++){
			
			double[] fArr = new double[2];
			PriorityQueue<Organism> a = c.get(keysC[i]);
			int size = a.size();
			
			for(int j = 0; j < size; j++){
				
				if(j == 0)
					fArr[0] = a.poll().getFitness();
				else if(j > size/2 - 1 && j < size/2 + 1)
					fArr[1] = a.poll().getFitness();
				else
					a.poll();
				
			}
			
			fitnessMetrics.put(keysC[i], fArr);
			
		}
		
		for(int i = 0; i < keysC.length; i++){
			
			if(extinctionMap.containsKey(keysC[i])){
				
				double[] a = fitnessMetrics.get(keysC[i]);
				double[] b = extinctionMap.get(keysC[i]);
				
				double penalty = 0;
				
				if(a[0] < b[0])
					penalty += 0.25;
				if(a[1] <= b[1])
					penalty += 0.75;
				
				b[2] += penalty;
				
				extinctionMap.put(keysC[i], b);
				
			}
			else{
				double[] adeeb = new double[3];
				adeeb[0] = fitnessMetrics.get(keysC[i])[0];
				adeeb[1] = fitnessMetrics.get(keysC[i])[1];
				adeeb[2] = 0;
				extinctionMap.put(keysC[i], adeeb);
			}
			
		}
		
		PriorityQueue<Organism> queue = new PriorityQueue<Organism>(fitnessComparator);

		for(int i = 0; i < prevGeneration.size(); i++)
			queue.add(prevGeneration.get(i));
		
		double avgFitness = 0;
		double avgFitnessU = 0;
		
		for(int i = 0; i < prevGeneration.size(); i++)
			if(i > prevGeneration.size()/2 - 1 && i < prevGeneration.size()/2 + 1){
			
				Organism med = queue.poll();
				avgFitness = med.getFitness();
				avgFitnessU = avgFitness*speciesSize.get(med.getSpecies());
				break;
				
			}
			else
				queue.poll();
		System.out.println("average fitness:" + Math.round(avgFitnessU));
		
		HashMap<Integer, Double> speciesFitnessThreshold = new HashMap<Integer, Double>();
		HashMap<Integer, Double> speciesNumOffspring = new HashMap<Integer, Double>();
		
		speciesSize.clear();
		
		for(int i = 0; i < prevGeneration.size(); i++)
			if(speciesSize.containsKey(prevGeneration.get(i).getSpecies()))
				speciesSize.put(prevGeneration.get(i).getSpecies(), speciesSize.get(prevGeneration.get(i).getSpecies()) + 1);
			else
				speciesSize.put(prevGeneration.get(i).getSpecies(), 1);
		
		Integer[] species = speciesSize.keySet().toArray(new Integer[0]);
		Organism[] organisms = prevGeneration.values().toArray(new Organism[0]);
		for(int i = 0; i < species.length; i++){
			
			int currSpecies = species[i];
			speciesNumOffspring.put(currSpecies, 0.0);
			PriorityQueue<Organism> fitnessQueue = new PriorityQueue<Organism>(fitnessComparator);
			
			for(int j = 0; j < organisms.length; j++)
				if(organisms[j].getSpecies() == currSpecies){
					
					fitnessQueue.add(organisms[j]);
					speciesNumOffspring.put(currSpecies, speciesNumOffspring.get(currSpecies) + organisms[j].getFitness()/avgFitness);
					
				}
			
			int size = fitnessQueue.size();
			for(int j = 0; j < size; j++)
				if(j <= (((double)size)*0.2 + 0.5) && j >= (((double)size)*0.2 - 0.5))
					speciesFitnessThreshold.put(currSpecies, fitnessQueue.poll().getFitness());
				else
					fitnessQueue.poll();
			
		}
		
		organisms = prevGeneration.values().toArray(new Organism[0]);
		for(int i = 0; i < organisms.length; i++)
			if(organisms[i].getFitness() < speciesFitnessThreshold.get(organisms[i].getSpecies()))
				prevGeneration.remove(organisms[i].getID());
		
		double nextPop = 0;
		
		for(int i = 0; i < organisms.length; i++)
			nextPop += Math.round(1.0*speciesNumOffspring.get(organisms[i].getSpecies())/(0.2*speciesSize.get(organisms[i].getSpecies())));
		
		//Cause extinction
		for(int i = 0; i < keysC.length; i++){
			
			if(extinctionMap.get(keysC[i])[2] > 15){
				
				organisms = prevGeneration.values().toArray(new Organism[0]);
				for(int j = 0; j < organisms.length; j++)
					if(organisms[j].getSpecies() == keysC[i])
						prevGeneration.remove(organisms[j].getID());
				
				Integer[] keysF = speciesNumOffspring.keySet().toArray(new Integer[0]);
				double popMove = 0;
				for(int j = 0; j < keysF.length; j++)
					if(keysC[i] == keysF[j]){
						popMove = speciesNumOffspring.get(keysF[j]);
						speciesNumOffspring.remove(keysF[j]);
						break;
					}
				
				for(int j = 0; j < keysF.length; j++){
					
					if(keysF[j] == keysC[i] && j != keysF.length - 1)
						j++;
					else if(keysF[j] == keysC[i] && j == keysF.length - 1)
						break;
					
					double adeeb = speciesNumOffspring.get(keysF[j]);
					
					speciesNumOffspring.put(keysF[j], adeeb + popMove*speciesNumOffspring.get(keysF[j])/nextPop);
					
				}
				
			}
			
		}
		
		double offspringFactor = 1.0;
		//fix population constraints
		if(nextPop > 500)
			offspringFactor = 0.8*50/(nextPop - 500) + 0.2;
		offspringFactor = Math.min(offspringFactor, 1.0);
		offspringFactor = Math.max(offspringFactor, 0.2);
		nextPop = 0;
		organisms = prevGeneration.values().toArray(new Organism[0]);
		for(int i = 0; i < organisms.length; i++){
			nextPop += Math.round(offspringFactor*speciesNumOffspring.get(organisms[i].getSpecies())/(0.2*speciesSize.get(organisms[i].getSpecies())));
		}
			
		while(nextPop > 1000){
			
			offspringFactor *= 0.7;
			nextPop = 0;
			for(int i = 0; i < organisms.length; i++)
				nextPop += Math.round(offspringFactor*speciesNumOffspring.get(organisms[i].getSpecies())/(0.2*speciesSize.get(organisms[i].getSpecies())));
			
		}
		
		System.out.println("predicted next population:" + nextPop);
		
		Random random = new Random();
		for(int i = 0; i < organisms.length; i++){
			
			Organism parent = organisms[i];
			int numChildren = (int) Math.round(offspringFactor*speciesNumOffspring.get(parent.getSpecies())/(0.2*speciesSize.get(parent.getSpecies())));
			
			for(int j = 0; j < numChildren; j++){
			
				int randID = random.nextInt(organisms.length);
				while(organisms[randID].getSpecies() != parent.getSpecies())
					randID = random.nextInt(organisms.length);
				Organism parent1 = organisms[randID];
				Organism child1 = makeChild(parent, parent1, true);
				child1.setID(numOrgo);
				generation.put(numOrgo, child1);
				numOrgo++;
			
			}
			
		}
		
		prevGeneration = new HashMap<Integer, Organism>(prevGenerationTemp);
		
		generations++;
		double[][] distanceMap = distanceMap(prevGeneration);
		speciate(distanceMap, prevGeneration);
		System.out.println("Generation:" + generations + ", Delta:" + delta + ", Organisms: " + generation.size() + ", Active Species: " + activeSpecies);
		fitness();
		
	}
	
	private double[][] distanceMap(HashMap<Integer, Organism> prevGeneration){
		
		double[][] distanceMap = new double[generation.size()][prevGeneration.size()];
		
		for(int i = 0; i < generation.size(); i++){
			
			for(int j = 0; j < prevGeneration.size(); j++){
				
				HashMap<Integer, Gene> genomeI = generation.get(i).cloneGenome();
				Integer[] keysI = genomeI.keySet().toArray(new Integer[0]);
				HashMap<Integer, Gene> genomeIEdges = new HashMap<Integer, Gene>(genomeI);
				for(int a = 0; a < keysI.length; a++)
					if(!genomeIEdges.get(keysI[a]).getIsEdge())
						genomeIEdges.remove(keysI[a]);
				keysI = genomeIEdges.keySet().toArray(new Integer[0]);
				
				HashMap<Integer, Gene> genomeJ = prevGeneration.get(j).cloneGenome();
				Integer[] keysJ = genomeJ.keySet().toArray(new Integer[0]);
				HashMap<Integer, Gene> genomeJEdges = new HashMap<Integer, Gene>(genomeJ);
				for(int b = 0; b < keysJ.length; b++)
					if(!genomeJEdges.get(keysJ[b]).getIsEdge())
						genomeJEdges.remove(keysJ[b]);
				keysJ = genomeJEdges.keySet().toArray(new Integer[0]);
				
				double N = (double)Math.max(keysI.length, keysJ.length);
				
				double Wn = 0;
				double n = 0;
				double D = 0;
				double E = 0;
				double maxI = 0;
				double maxJ = 0;
				
				for(int a = 0; a < genomeIEdges.size(); a++)
					if(genomeIEdges.get(keysI[a]).getIsEdge())
						maxI = (double) Math.max(keysI[a], maxI);
				for(int b = 0; b < genomeJEdges.size(); b++)
					if(genomeJEdges.get(keysJ[b]).getIsEdge())
						maxJ = (double) Math.max(keysJ[b], maxJ);
				
				for(int a = 0; a < keysI.length; a++){
					
					for(int b = 0; b < keysJ.length; b++){
						
						if(keysI[a] == keysJ[b]){
							
							n++;
							Wn += Math.abs(genomeI.get(keysI[a]).getWeight() - genomeJ.get(keysJ[b]).getWeight());
							
						}
						else if(genomeI.get(keysI[a]).getIn() == genomeJ.get(keysJ[b]).getIn() && genomeI.get(keysI[a]).getOut() == genomeJ.get(keysJ[b]).getOut()){
							
							n++;
							Wn += Math.abs(genomeI.get(keysI[a]).getWeight() - genomeJ.get(keysJ[b]).getWeight());
							
						}
						else if(keysI[a] > maxJ && b == 0)
							E++;
						else if(keysJ[b] > maxI && a == 0)
							E++;
						
					}
					
				}
				
				D = Math.min(keysI.length, keysJ.length) - n;
				
				double W;
				
				if(n == 0)
					W = Double.MAX_VALUE;
				else
					W = Wn/n;
				
				distanceMap[i][j] = E/N + D/N + 0.2*W;
				
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
	
			if(i == (int)(size*0.2))
				delta = queue2.poll();
			else
				queue2.poll();
			
		}
		
		delta = Math.min(delta, 0.2 + 0.3/generations);
		
		return distanceMap;
		
	}
	
	private double[][] distanceMap(){
		
		double[][] distanceMap = new double[generation.size()][generation.size()];
		
		for(int i = 0; i < generation.size(); i++){
			
			for(int j = 0; j < generation.size(); j++){
				
				if(i != j && distanceMap[j][i] == 0){
					
					HashMap<Integer, Gene> genomeI = generation.get(i).cloneGenome();
					Integer[] keysI = genomeI.keySet().toArray(new Integer[0]);
					HashMap<Integer, Gene> genomeIEdges = new HashMap<Integer, Gene>(genomeI);
					for(int a = 0; a < keysI.length; a++)
						if(!genomeIEdges.get(keysI[a]).getIsEdge())
							genomeIEdges.remove(keysI[a]);
					keysI = genomeIEdges.keySet().toArray(new Integer[0]);
					
					HashMap<Integer, Gene> genomeJ = generation.get(j).cloneGenome();
					Integer[] keysJ = genomeJ.keySet().toArray(new Integer[0]);
					HashMap<Integer, Gene> genomeJEdges = new HashMap<Integer, Gene>(genomeJ);
					for(int b = 0; b < keysJ.length; b++)
						if(!genomeJEdges.get(keysJ[b]).getIsEdge())
							genomeJEdges.remove(keysJ[b]);
					keysJ = genomeJEdges.keySet().toArray(new Integer[0]);
					
					int N = Math.max(keysI.length, keysJ.length);
					
					double Wn = 0;
					double n = 0;
					double D = 0;
					double E = 0;
					double maxI = 0;
					double maxJ = 0;
					
					for(int a = 0; a < genomeIEdges.size(); a++)
						if(genomeI.get(keysI[a]).getIsEdge())
							maxI = (double) Math.max(keysI[a], maxI);
					for(int b = 0; b < genomeJEdges.size(); b++)
						if(genomeJ.get(keysJ[b]).getIsEdge())
							maxJ = (double) Math.max(keysJ[b], maxJ);
					
					for(int a = 0; a < keysI.length; a++){
						
						for(int b = 0; b < keysJ.length; b++){
							
							if(keysI[a] == keysJ[b]){
								
								n++;
								Wn += Math.abs(genomeI.get(keysI[a]).getWeight() - genomeJ.get(keysJ[b]).getWeight());
								
							}
							else if(genomeI.get(keysI[a]).getIn() == genomeJ.get(keysJ[b]).getIn() && genomeI.get(keysI[a]).getOut() == genomeJ.get(keysJ[b]).getOut()){
								
								n++;
								Wn += Math.abs(genomeI.get(keysI[a]).getWeight() - genomeJ.get(keysJ[b]).getWeight());
								
							}
							else if(keysI[a] > maxJ && b == 0)
								E++;
							else if(keysJ[b] > maxI && a == 0)
								E++;
							
						}
						
					}
					
					D = Math.min(keysI.length, keysJ.length) - n;
					
					double W;
					
					if(n == 0)
						W = Double.MAX_VALUE;
					else
						W = Wn/n;
					
					distanceMap[i][j] = E/N + D/N + W;
					distanceMap[j][i] = distanceMap[i][j];
					if(E/N > 1 || D/N > 1)
						System.out.println("E:" + E + " D:" + D + " W:" + W + " N:" + N);
					
				}
				
			}
			
		}
		
		return distanceMap;
		
	}
	
	private void speciate(double[][] distanceMap, HashMap<Integer, Organism> prevGeneration){
		
		activeSpecies = 0;
		HashMap<Integer, Integer> prevSpeciesSize = new HashMap<Integer, Integer>(speciesSize);
		int maxPrev = Collections.max(prevSpeciesSize.keySet());
		nextSpecies = Math.max(maxPrev + 1, nextSpecies);
		double[][] distanceMapAlt = distanceMap();
		speciesSize.clear();
		
		for(int i = 0; i < generation.size(); i++){

			HashMap<Integer, Double> avgSpeciesDist = new HashMap<Integer, Double>();
			for(int j = 0; j < prevGeneration.size(); j++)
				if(avgSpeciesDist.containsKey(prevGeneration.get(j).getSpecies()))
					avgSpeciesDist.put(prevGeneration.get(j).getSpecies(), avgSpeciesDist.get(prevGeneration.get(j).getSpecies()) + distanceMap[i][j]*prevGeneration.get(j).getFitness()/prevSpeciesSize.get(prevGeneration.get(j).getSpecies()));
				else
					avgSpeciesDist.put(prevGeneration.get(j).getSpecies(), distanceMap[i][j]*prevGeneration.get(j).getFitness()/prevSpeciesSize.get(prevGeneration.get(j).getSpecies()));
			Integer[] keys = avgSpeciesDist.keySet().toArray(new Integer[0]);
			for(int j = 0; j < keys.length; j++)
				if(prevSpeciesSize.get(keys[j]) > 0)
					avgSpeciesDist.put(keys[j], avgSpeciesDist.get(keys[j])/(double)prevSpeciesSize.get(keys[j]));
			keys = shuffleArray(keys);
			for(int j = 0; j < keys.length; j++){
				
				if(avgSpeciesDist.get(keys[j]) < delta){
					
					generation.get(i).setSpecies(keys[j]);
					if(speciesSize.containsKey(keys[j]))
						speciesSize.put(keys[j], speciesSize.get(keys[j]) + 1);
					else{
						
						speciesSize.put(keys[j], 1);
						activeSpecies++;
						
					}
					break;
					
				}
				else if(j == keys.length - 1 && i == 0){
					
					generation.get(i).setSpecies(nextSpecies);
					speciesSize.put(nextSpecies, 1);
					nextSpecies++;
					activeSpecies++;
					
				}
				else if(j == keys.length - 1 && i != 0){
					
					avgSpeciesDist = new HashMap<Integer, Double>();
					for(int k = 0; k < i; k++){
						if(avgSpeciesDist.containsKey(generation.get(k).getSpecies()))
							avgSpeciesDist.put(generation.get(k).getSpecies(), avgSpeciesDist.get(generation.get(k).getSpecies()) + distanceMapAlt[i][k]);
						else
							avgSpeciesDist.put(generation.get(k).getSpecies(), distanceMapAlt[i][k]);
					
						if(generation.get(k).getSpecies() == -1)
							System.out.println(k + " " + i);
					}
					keys = avgSpeciesDist.keySet().toArray(new Integer[0]);
					for(int k = 0; k < keys.length; k++)
						if(speciesSize.get(keys[k]) == null)
							System.out.println(keys[k] + " " + keys.length + " " + speciesSize.size());
					for(int k = 0; k < keys.length; k++)
						if(speciesSize.get(keys[k]) > 0)
							avgSpeciesDist.put(keys[k], avgSpeciesDist.get(keys[k])/(double)speciesSize.get(keys[k]));
					keys = shuffleArray(keys);
					for(int k = 0; k < keys.length; k++){
						
						if(avgSpeciesDist.get(keys[k]) < delta && keys[k] > maxPrev){
							
							generation.get(i).setSpecies(keys[k]);
							speciesSize.put(keys[k], speciesSize.get(keys[k]) + 1);
							break;
							
						}
						else if(k == keys.length - 1){
							
							generation.get(i).setSpecies(nextSpecies);
							speciesSize.put(nextSpecies, 1);
							nextSpecies++;
							activeSpecies++;
							
						}
						
					}
					
				}
				
			}
			
		}

		speciesSize.clear();
		for(int i = 0; i < generation.size(); i++)
			if(speciesSize.containsKey(generation.get(i).getSpecies()))
				speciesSize.put(generation.get(i).getSpecies(), speciesSize.get(generation.get(i).getSpecies()) + 1);
			else
				speciesSize.put(generation.get(i).getSpecies(), 1);
			
		/*Random random = new Random();
		Integer[] species = speciesSize.keySet().toArray(new Integer[0]);
		for(int k = 0; k < speciesSize.size(); k++){

			//if there are an odd number of organisms
			//randomly breed
			int i = species[k];
			if(speciesSize.get(i) % 2 == 1 && speciesSize.get(i) > 0){
				
				int randID1 = random.nextInt(generation.size());
				while(generation.get(randID1).getSpecies() != i)
					randID1 = random.nextInt(generation.size());
					
				int randID2;
				if(speciesSize.get(i) > 1){
					randID2 = random.nextInt(generation.size());
					while(generation.get(randID2).getSpecies() != i)
						randID2 = random.nextInt(generation.size());
						
				}
				else
					randID2 = randID1;
				Organism parent1 = generation.get(randID1);
				Organism parent2 = generation.get(randID2);
				Organism child = makeChild(parent1, parent2, false);
				child.setID(generation.size());
				child.setSpecies(i);
				generation.put(generation.size(), child);
				
				speciesSize.put(i, speciesSize.get(i) + 1);
				
			}
			
		}

		for(int k = 0; k < speciesSize.size(); k++){
			
			int i = species[k];
			while(speciesSize.get(i) < 4 && speciesSize.get(i) > 0){
				
				int randID1 = random.nextInt(generation.size());
				while(generation.get(randID1).getSpecies() != i)
					randID1 = random.nextInt(generation.size());
				int randID2;
				if(speciesSize.get(i) > 1){
					randID2 = random.nextInt(generation.size());
					while(generation.get(randID2).getSpecies() != i)
						randID2 = random.nextInt(generation.size());
				}
				else
					randID2 = randID1;
				Organism parent1 = generation.get(randID1);
				Organism parent2 = generation.get(randID2);
				Organism child = makeChild(parent1, parent2, false);
				child.setID(generation.size());
				child.setSpecies(i);
				generation.put(generation.size(), child);
				
				speciesSize.put(i, speciesSize.get(i) + 1);
				
			}
			
		}*/
		
	}

	private Organism makeChild(Organism parent1, Organism parent2, boolean canMutate){
		
		HashMap<Integer, Gene> genome1 = parent1.cloneGenome();
		Gene[] genome1values = genome1.values().toArray(new Gene[0]);
		HashMap<Integer, Gene> genome2 = parent2.cloneGenome();
		Gene[] genome2values = genome2.values().toArray(new Gene[0]);
		
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
		
		double willMutate = random.nextDouble();
		if(canMutate && willMutate < 0.5)
			genome = mutate(genome);
		
		return new Organism(genome);
		
	}
	
	
	private HashMap<Integer, Gene> mutate(HashMap<Integer, Gene> genome){
		
		Random random = new Random();
		double numMutations = Math.abs(random.nextGaussian())*1.5 + 1.0;
		for(int i = 0; i < numMutations; i++){
			
			int mutation = random.nextInt(4);
			if(mutation <= 1)
				genome = addNode(genome);
			if(mutation == 2)
				genome = addEdge(genome);
			if(mutation == 3)
			genome = adjustWeight(genome);
			
		}
		
		return genome;
		
	}

	private HashMap<Integer, Gene> addNode(HashMap<Integer, Gene> genome){
		
		int id = 0;
		
		Gene[] genomeArr = genome.values().toArray(new Gene[0]);
		ArrayList<Gene> edges = new ArrayList<Gene>();
		for(int i = 0; i < genomeArr.length; i++)
			if(genomeArr[i].getIsEdge() && genomeArr[i].getEnabled()){
				
				if(!edges.contains(genomeArr[i]))
					edges.add(genomeArr[i]);
				
			}
		
		if(edges.size() == 0){
			for(int i = 0; i < genomeArr.length; i++)
				if(genomeArr[i].getIsEdge()){
					genomeArr[i].enable();
					genome.put(genomeArr[i].getInnov(), genomeArr[i]);
				}
			return genome;
		}
		
		Random random = new Random();
		Gene edge = edges.get(random.nextInt(edges.size()));
		genome.remove(edge.getInnov());
		edge.disable();
		genome.put(edge.getInnov(), edge);
		int rand1 = edge.getIn();
		int rand2 = edge.getOut();
		
		Gene node1 = new Gene(-1, false, -1, -1, -1, -1, -1);
		Gene node2 = new Gene(-1, false, -1, -1, -1, -1, -1);
		
		for(int i = 0; i < genomeArr.length; i++)
			if(!genomeArr[i].getIsEdge()){
				
				if(genomeArr[i].getID() == rand1)
					node1 = genomeArr[i];
				else if(genomeArr[i].getID() == rand2)
					node2 = genomeArr[i];
				
			}
		
		Gene[] geneListArr = geneList.values().toArray(new Gene[0]);
		for(int i = 0; i < geneListArr.length; i++)
			id = Math.max(id, geneListArr[i].getID());
		id++;
		
		double gaussian = random.nextGaussian();
		double tierFactor = 2.0;
		
		if(gaussian >= 1.0)
			tierFactor = 3.0;
		else if(gaussian <= -1.0)
			tierFactor = 1.0;
		
		Gene node = new Gene(innovTotal, false, -1, -1, -1, id, tierFactor*(node1.getTier() + node2.getTier())/4);
		Gene edge1 = new Gene(innovTotal + 1, true, node.getID(), rand2, 1.0, -1, -1);
		Gene edge2 = new Gene(innovTotal + 2, true, rand1, node.getID(), 1.0, -1, -1);
		
		top:
		for(int i = 0; i < geneListArr.length; i++){
			
			if(geneListArr[i].getIsEdge() && geneListArr[i].getOut() == rand2){
				
				int nodeID = geneListArr[i].getIn();
				
				for(int j = 0; j < geneListArr.length; j++){
					
					if(j == i && i != geneListArr.length - 1)
						j++;
					
					if(geneListArr[j].getIsEdge() && geneListArr[j].getOut() == nodeID && geneListArr[j].getIn() == rand1){
						
						for(int k = 0; k < geneListArr.length; k++){
							
							if(k == i && i != geneListArr.length - 1)
								k++;
							if(k == j && j != geneListArr.length - 1)
								k++;
							
							if(geneListArr[k].getTier() == node.getTier() && geneListArr[k].getID() == nodeID){
								
								double doOrDie = random.nextDouble();
								
								if(doOrDie > 0.25){
									
									node = geneListArr[k].clone();
									edge1 = geneListArr[i].clone();
									edge2 = geneListArr[j].clone();
									
								}
								
								break top;
								
							}
							
						}
						
					}
					
				}
				
			}
			
		}
		
		if(node.getInnov() >= innovTotal){
			
			geneList.put(node.getInnov(), node);
			innovTotal++;
			
			if(edge1.getInnov() >= innovTotal){
				
				geneList.put(edge1.getInnov(), edge1);
				innovTotal++;
				
				if(edge2.getInnov() >= innovTotal){
					
					geneList.put(edge2.getInnov(), edge2);
					innovTotal++;
					
				}
				
			}
			
		}
		
		node = node.clone();
		edge1 = edge1.clone();
		edge1.setWeight(edge.getWeight());
		edge2 = edge2.clone();
		edge2.setWeight(edge.getWeight());
		
		genome.put(node.getInnov(), node);
		genome.put(edge1.getInnov(), edge1);
		genome.put(edge2.getInnov(), edge2);
		
		return genome;
		
	}

	private HashMap<Integer, Gene> addEdge(HashMap<Integer, Gene> genome){
		
		ArrayList<Gene> nodes = new ArrayList<Gene>();
		ArrayList<Gene> edges = new ArrayList<Gene>();
		Gene[] genomeArr = genome.values().toArray(new Gene[0]);
		for(int i = 0; i < genomeArr.length; i++)
			if(!genomeArr[i].getIsEdge() && !nodes.contains(genomeArr[i]))
				nodes.add(genomeArr[i]);
			else if(!edges.contains(genomeArr[i]))
				edges.add(genomeArr[i]);
		
		Random random = new Random();
		
		int rand1 = random.nextInt(nodes.size());
		int rand2 = random.nextInt(nodes.size());
		int edge = -1;
		
		for(int i = 0; i < edges.size(); i++)
			if(edges.get(i).getIn() == nodes.get(rand1).getID() && edges.get(i).getOut() == nodes.get(rand2).getID())
				edge = i;
		
		boolean full = true;
		
		for(int i = 0; i < nodes.size(); i++){
			
			for(int j = 0; j < nodes.size(); j++){
				
				if(nodes.get(i).getTier() > nodes.get(j).getTier()){
					
					boolean exists = false;
					
					for(int k = 0; k < edges.size(); k++){
						
						if(edges.get(k).getIn() == nodes.get(i).getID() && edges.get(k).getOut() == nodes.get(j).getID())
							exists = true;
						
					}
					
					full = full && exists;
					
				}
				
			}
			
		}
		
		if(full)
			return genome;
		
		while(nodes.get(rand1).getTier() <= nodes.get(rand2).getTier() || edge >= 0){
			
			rand1 = random.nextInt(nodes.size());
			rand2 = random.nextInt(nodes.size());
			edge = -1;
			
			for(int i = 0; i < edges.size(); i++)
				if(edges.get(i).getIn() == nodes.get(rand1).getID() && edges.get(i).getOut() == nodes.get(rand2).getID())
					edge = i;
			
		}
		
		Gene gene = new Gene(innovTotal, true, nodes.get(rand1).getID(), nodes.get(rand2).getID(), 1.0, -1, -1);
		
		Gene[] geneListArr = geneList.values().toArray(new Gene[0]);
		for(int i = 0; i < geneListArr.length; i++)
			if(geneListArr[i].getIn() == gene.getIn() && geneListArr[i].getOut() == gene.getOut())
				gene = geneListArr[i].clone();
		
		if(gene.getInnov() >= innovTotal){
			
			geneList.put(gene.getInnov(), gene);
			gene = gene.clone();
			innovTotal++;
			
		}
		
		gene = gene.clone();
		gene.setWeight(2*random.nextDouble() - 1);
		genome.put(gene.getInnov(), gene);
		
		return genome;
		
	}
	
	private HashMap<Integer, Gene> adjustWeight(HashMap<Integer, Gene> genome){
		
		ArrayList<Gene> edges = new ArrayList<Gene>();
		Gene[] genomeArr = genome.values().toArray(new Gene[0]);
		for(int i = 0; i < genomeArr.length; i++)
			if(genomeArr[i].getIsEdge() && !edges.contains(genomeArr[i]))
				edges.add(genomeArr[i]);
		
		Random random = new Random();
		Gene edge = edges.get(random.nextInt(edges.size()));
		genome.remove(edge.getInnov());
		
		double chance = random.nextDouble();
		
		double weight = edge.getWeight() + 0.25*random.nextGaussian();
		if(weight > 2)
			weight = 2;
		else if(weight < -2)
			weight = -2;
		
		if(chance >= 0.2)
			edge.setWeight(weight);
		else
			edge.setWeight(2.0*random.nextDouble() - 1.0);
		
		genome.put(edge.getInnov(), edge);
		
		return genome;
		
	}
	
	private void fitness(){
		
		double maxfit = 0;
		int maxfitspec = 0;
		double minfit = Double.MAX_VALUE;
		int minfitspec = 0;
		int n = 30;

		for(int i = 0; i < generation.size(); i++){
			
			Brain brain = genoToPheno(generation.get(i).cloneGenome());
			
			for(int j = 0; j < n; j++){
				
				int score = 0;
				theGame:
				while(game.canWin){
					
					int[][] boardB4 = new int[game.getBoard().length][game.getBoard()[0].length];
					for(int a = 0; a < game.getBoard().length; a++)
						for(int b = 0; b < game.getBoard()[0].length; b++)
							boardB4[a][b] = game.getBoard()[a][b];
					
					score = game.score;
					double[] moveset = brain.brainMove(game.getBoard());
					HashMap<Double, Integer> thing = new HashMap<Double, Integer>();
					int[] moveorder = new int[4];
					
					for(int a = 0; a < 4; a++)
						thing.put(moveset[a], a);
					
					Arrays.sort(moveset);
					
					for(int a = 0; a < 4; a++)
						moveorder[a] = thing.get(moveset[a]);
					
					boolean equals = true;
					for(int a = 0; a < game.getBoard().length; a++)
						for(int b = 0; b < game.getBoard()[a].length; b++)
							equals = equals && (boardB4[a][b] == game.getBoard()[a][b]);
		
					//If it doesn't make a move that changes anything
					//while(boardB4.equals(game.getBoard()))
					//	game.move(random.nextInt(4));
					int c = 0;
					
					while(equals){
						
						if(c == 4)
							break theGame;
						
						game.move(moveorder[c]);
						
						c++;
						
						equals = true;
						for(int a = 0; a < game.getBoard().length; a++)
							for(int b = 0; b < game.getBoard()[a].length; b++)
								equals = equals && (boardB4[a][b] == game.getBoard()[a][b]);
					}
					
				}
				score += generation.get(i).getFitness();
				generation.get(i).setFitness((double)score);
				game.resetBoard();
				
			}
			
			generation.get(i).setFitness(generation.get(i).getFitness()/n);
			if(generation.get(i).getFitness() > maxfit){
				
				maxfit = generation.get(i).getFitness();
				maxfitspec = i;
				
			}
			
			if(generation.get(i).getFitness() < minfit){
				
				minfit = generation.get(i).getFitness();
				minfitspec = i;
				
			}
			generation.get(i).setFitness(generation.get(i).getFitness()/((double)speciesSize.get(generation.get(i).getSpecies())));
			
		}
		
		maxFitness = maxfit;
		
		System.out.println("Organism " + maxfitspec + ", Species " + generation.get(maxfitspec).getSpecies() + "\nfitness:" + Math.round(maxFitness));
		
		double[][] fit2 = genoToPheno(generation.get(maxfitspec).cloneGenome()).getGraph();
		System.out.println("nodes:" + fit2.length);
		
		if(maxFitness > 2000){
		
			System.out.println("topology:");
			double[][] fit = genoToPheno(generation.get(maxfitspec).cloneGenome()).getGraph();
			
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
		
	}
	
	
	private Brain genoToPheno(HashMap<Integer, Gene> genome){
		
 		Gene[] genes = genome.values().toArray(new Gene[0]);
		
		HashMap<Integer, Double> nodesH = new HashMap<Integer, Double>();
		Gene[] genomeKeys = genome.values().toArray(new Gene[0]);
		
		for(int i = 0; i < 4; i++)
			nodesH.put(i, Double.MAX_VALUE/2);
		for(int i = 4; i < 20; i++)
			nodesH.put(i, 0.0);
		for(int i = 0; i < genomeKeys.length; i++){
			
			if(!genomeKeys[i].getIsEdge()){
				if(!nodesH.containsKey(genomeKeys[i].getID()))
					nodesH.put(genomeKeys[i].getID(), genomeKeys[i].getTier());
					
			}
			
		}

		PriorityQueue<Double> nodesL = new PriorityQueue<Double>();
		Integer[] vertices = nodesH.keySet().toArray(new Integer[0]);
		Integer[] verticesClone = new Integer[vertices.length];
		for(int i = 0; i < vertices.length; i++)
			verticesClone[i] = vertices[i];
		
		for(int i = 0; i < vertices.length; i++)
			nodesL.add(nodesH.get(vertices[i]));
		
		for(int i = 0; i < vertices.length; i++){
			
			double tier = nodesL.poll();
			
			for(int j = 0; j < verticesClone.length; j++){
				
				if(tier == nodesH.get(verticesClone[j]))
					vertices[i] = verticesClone[j];
				
			}
			
		}
		
		double[][] graph = new double[nodesH.size()][nodesH.size()];
		
		for(int i = 0; i < genes.length; i++){
			
			if(genes[i].getEnabled() && genes[i].getIsEdge()){
				
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
	
	
	public Integer[] shuffleArray(Integer[] ar){
		
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
