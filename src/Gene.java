public class Gene {
	
	private int innovNum;
	private boolean isEdge;
	private int in;
	private int out;
	private double weight;
	private int id;
	private double tier;
	private boolean enabled;
	
	public Gene(int innovNum, boolean isEdge, int in, int out, double weight, int id, double tier){
		
		this.innovNum = innovNum;
		this.isEdge = isEdge;
		if(isEdge){
			
			this.in = in;
			this.out = out;
			this.weight = weight;
			this.id = -1;
			this.tier = -1;
			
		}
		else{
			
			this.in = -1;
			this.out = -1;
			this.weight = 0;
			this.id = id;
			this.tier = tier;
			
		}
		this.enabled = true;
		
	}
	
	public int getInnov(){
		
		return innovNum;
		
	}
	
	public boolean getIsEdge(){
		
		return isEdge;
		
	}
	
	public int getIn(){
		
		return in;
		
	}
	
	public int getOut(){
		
		return out;
		
	}
	
	public double getWeight(){
		
		return weight;
		
	}
	
	public int getID(){
		
		return id;
		
	}
	
	public double getTier(){
		
		return tier;
		
	}
	
	public boolean getEnabled(){
		
		return enabled;
		
	}
	
	public void disable(){
		
		enabled = false;
		
	}
	
	public void enable(){
		
		enabled = true;
		
	}

	public void setWeight(double weight){
		
		this.weight = weight;
		
	}
	
	public void setInnov(int innov){
		
		this.innovNum = innov;
		
	}
	
	public Gene clone(){
		
		Gene gene = new Gene(innovNum, isEdge, in, out, weight, id, tier);
		if(!enabled)
			gene.disable();
		
		return gene;
		
	}
	
}