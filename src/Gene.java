public class Gene {
	
	private int innovNum;
	private int in;
	private int out;
	private double weight;
	private boolean enabled;
	
	public Gene(int in, int out, double weight, int innovNum){
		
		this.in = in;
		this.out = out;
		this.weight = weight;
		this.innovNum = innovNum;
		enabled = true;
		
	}
	
	public boolean getEnabled(){
		
		return enabled;
		
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
	
	public int getInnov(){
		
		return innovNum;
		
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
		
		Gene gene = new Gene(in, out, weight, innovNum);
		if(!enabled)
			gene.disable();
		
		return gene;
		
	}
	
}