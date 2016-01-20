import java.util.*;

public class Organism{
	
	private HashMap<Integer, Gene> genome;
	private int id;
	private double fitness;
	private int species;
	
	public Organism(){
		
		id = -1;
		fitness = -1;
		species = -1;
		
		genome = new HashMap<Integer, Gene>();
		
	}
	
	public Organism(HashMap<Integer, Gene> genome){
		
		this.id = -1;
		fitness = -1;
		species = -1;
		
		this.genome = new HashMap<Integer, Gene>(genome);
		
	}
	
	public Organism(int id){
		
		this.id = id;
		fitness = -1;
		species = -1;
		
		genome = new HashMap<Integer, Gene>();
		
	}
	
	public Organism(int id, HashMap<Integer, Gene> genome){
		
		this.id = id;
		fitness = -1;
		species = -1;
		
		this.genome = new HashMap<Integer, Gene>(genome);
		
	}
	
	public int getID(){
		
		return id;
		
	}
	
	public double getFitness(){
		
		return fitness;
		
	}
	
	public int getSpecies(){
		
		return species;
		
	}
	
	public HashMap<Integer, Gene> cloneGenome(){
		
		HashMap<Integer, Gene> clone = new HashMap<Integer, Gene>(genome);
		return clone;
		
	}
	
	public Organism clone(){
		
		Organism clone = new Organism();
		clone.setGenome(this.cloneGenome());
		clone.setFitness(fitness);
		clone.setSpecies(species);
		
		return clone;
		
	}
	
	public void addGene(Gene gene){
		
		if(genome.containsKey(gene.getInnov())){
			
			genome.remove(gene.getInnov());
			genome.put(gene.getInnov(), gene.clone());
			
		}
		
	}
	
	public void removeGene(int id){
		
		if(genome.containsKey(id))
			genome.remove(id);
		
	}

	public void setGenome(HashMap<Integer, Gene> genome){
		
		this.genome = genome;
		
	}
	
	public void setID(int id){
		
		this.id = id;
		
	}
	
	public void setFitness(double fitness){
		
		this.fitness = fitness;
		
	}
	
	public void setSpecies(int species){
		
		this.species = species;
		
	}
	
}
