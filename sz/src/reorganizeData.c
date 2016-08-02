#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "sz.h"

void reorganizeData_float(float **newData, float* oriData, int r5, int r4, int r3, int r2, int r1)
{
	if(r2==0)
	{
		*newData = oriData;
	}
	else
	{
		int index = 0;
		int dataLength = computeDataLength(r5,r4,r3,r2,r1);		
		*newData = (float*)malloc(sizeof(float)*dataLength);
		float *p = *newData;
		if(r3==0)
		{
			int a2,a1,i2,i1,n2,n1,x1;
			for(a2=0;a2<r2;a2+=reOrgSize)
			{
				n2 = a2+reOrgSize;
				if(n2 > r2)
					n2 = r2;
				for(a1=0;a1<r1;a1+=reOrgSize)
				{
					n1 = a1+reOrgSize;
					if(n1 > r1)
						n1 = r1;
					for(i2=a2;i2<n2;i2++)
					{
						x1 = i2*r1;
						memcpy(p, &(oriData[x1+a1]), (n1-a1)*sizeof(float));
						p+=(n1-a1);
					}			
				}
			}
		}
		else if(r4==0)
		{
			int a3,a2,a1,i3,i2,i1,n3,n2,n1,x2,x1;
			for(a3=0;a3<r3;a3+=reOrgSize)
			{
				n3 = a3+reOrgSize;
				if(n3 > r3)
					n3 = r3;
				for(a2=0;a2<r2;a2+=reOrgSize)
				{
					n2 = a2+reOrgSize;
					if(n2 > r2)
						n2 = r2;
					for(a1=0;a1<r1;a1+=reOrgSize)
					{
						n1 = a1+reOrgSize;
						if(n1 > r1)
							n1 = r1;
						for(i3=a3;i3<n3;i3++)
						{
							x2 = i3*r2*r1;
							for(i2=a2;i2<n2;i2++)
							{
								x1 = i2*r1;
								//for(i1=a1;i1<n1;i1++)
								//	*(p++) = oriData[x2+x1+i1];
								memcpy(p, &(oriData[x2+x1+a1]), (n1-a1)*sizeof(float));
								p+=(n1-a1);
							}							
						}	
					}
				}				
			}	
		}
		else if(r5==0)
		{
			int a4,a3,a2,a1,i4,i3,i2,i1,n4,n3,n2,n1,x3,x2,x1,y2;
			for(a4=0;a4<r4;a4+=reOrgSize)
			{
				n4 = a4+reOrgSize;
				for(a3=0;a3<r3;a3+=reOrgSize)
				{
					n3 = a3+reOrgSize;
					if(n3 > r3)
						n3 = r3;
					for(a2=0;a2<r2;a2+=reOrgSize)
					{
						n2 = a2+reOrgSize;
						if(n2 > r2)
							n2 = r2;
						for(a1=0;a1<r1;a1+=reOrgSize)
						{
							n1 = a1+reOrgSize;
							if(n1 > r1)
								n1 = r1;
							for(i4=a4;i4<n4;i4++)
							{
								y2 = r2*r1;
								x3 = i4*r3*y2;
								for(i3=a3;i3<n3;i3++)
								{
									x2 = i3*y2;
									for(i2=a2;i2<n2;i2++)
									{
										x1 = i2*r1;
										memcpy(p, &(oriData[x3+x2+x1+a1]), (n1-a1)*sizeof(float));
										p+=(n1-a1);
									}							
								}								
							}
	
						}
					}				
				}
			}			
		}
		else
		{
			int a5,a4,a3,a2,a1,i5,i4,i3,i2,i1,n5,n4,n3,n2,n1,x4,x3,x2,x1,y3,y2;
			for(a5=0;a5<r5;a5+=reOrgSize)
			{
				n5 = a5+reOrgSize;
				for(a4=0;a4<r4;a4+=reOrgSize)
				{
					n4 = a4+reOrgSize;
					for(a3=0;a3<r3;a3+=reOrgSize)
					{
						n3 = a3+reOrgSize;
						if(n3 > r3)
							n3 = r3;
						for(a2=0;a2<r2;a2+=reOrgSize)
						{
							n2 = a2+reOrgSize;
							if(n2 > r2)
								n2 = r2;
							for(a1=0;a1<r1;a1+=reOrgSize)
							{
								n1 = a1+reOrgSize;
								if(n1 > r1)
									n1 = r1;
								for(i5=a5;i5<n5;i5++)
								{
									y2 = r2*r1;
									y3 = r3*y2;
									x4 = i5*r4*y3;
									for(i4=a4;i4<n4;i4++)
									{
										x3 = (a4+i4)*y3;
										for(i3=a3;i3<n3;i3++)
										{
											x2 = i3*y2;
											for(i2=a2;i2<n2;i2++)
											{
												x1 = i2*r1;
												memcpy(p, &(oriData[x4+x3+x2+x1+a1]), (n1-a1)*sizeof(float));
												p+=(n1-a1);
											}							
										}								
									}
								}
							}
						}				
					}
				}				
			}
		}
	}

}

void reorganizeData_double(double **newData, double* oriData, int r5, int r4, int r3, int r2, int r1)
{
	if(r2==0)
	{
		*newData = oriData;
	}
	else
	{
		int index = 0;
		int dataLength = computeDataLength(r5,r4,r3,r2,r1);		
		*newData = (double*)malloc(sizeof(double)*dataLength);
		double *p = *newData;
		if(r3==0)
		{
			int a2,a1,i2,i1,n2,n1,x1;
			for(a2=0;a2<r2;a2+=reOrgSize)
			{
				n2 = a2+reOrgSize;
				if(n2 > r2)
					n2 = r2;
				for(a1=0;a1<r1;a1+=reOrgSize)
				{
					n1 = a1+reOrgSize;
					if(n1 > r1)
						n1 = r1;
					for(i2=a2;i2<n2;i2++)
					{
						x1 = i2*r1;
						memcpy(p, &(oriData[x1+a1]), (n1-a1)*sizeof(double));
						p+=(n1-a1);
					}			
				}
			}
		}
		else if(r4==0)
		{
			int a3,a2,a1,i3,i2,i1,n3,n2,n1,x2,x1;
			for(a3=0;a3<r3;a3+=reOrgSize)
			{
				n3 = a3+reOrgSize;
				if(n3 > r3)
					n3 = r3;
				for(a2=0;a2<r2;a2+=reOrgSize)
				{
					n2 = a2+reOrgSize;
					if(n2 > r2)
						n2 = r2;
					for(a1=0;a1<r1;a1+=reOrgSize)
					{
						n1 = a1+reOrgSize;
						if(n1 > r1)
							n1 = r1;
						for(i3=a3;i3<n3;i3++)
						{
							x2 = i3*r2*r1;
							for(i2=a2;i2<n2;i2++)
							{
								x1 = i2*r1;
								//for(i1=a1;i1<n1;i1++)
								//	*(p++) = oriData[x2+x1+i1];
								memcpy(p, &(oriData[x2+x1+a1]), (n1-a1)*sizeof(double));
								p+=(n1-a1);
							}							
						}	
					}
				}				
			}	
		}
		else if(r5==0)
		{
			int a4,a3,a2,a1,i4,i3,i2,i1,n4,n3,n2,n1,x3,x2,x1,y2;
			for(a4=0;a4<r4;a4+=reOrgSize)
			{
				n4 = a4+reOrgSize;
				for(a3=0;a3<r3;a3+=reOrgSize)
				{
					n3 = a3+reOrgSize;
					if(n3 > r3)
						n3 = r3;
					for(a2=0;a2<r2;a2+=reOrgSize)
					{
						n2 = a2+reOrgSize;
						if(n2 > r2)
							n2 = r2;
						for(a1=0;a1<r1;a1+=reOrgSize)
						{
							n1 = a1+reOrgSize;
							if(n1 > r1)
								n1 = r1;
							for(i4=a4;i4<n4;i4++)
							{
								y2 = r2*r1;
								x3 = i4*r3*y2;
								for(i3=a3;i3<n3;i3++)
								{
									x2 = i3*y2;
									for(i2=a2;i2<n2;i2++)
									{
										x1 = i2*r1;
										memcpy(p, &(oriData[x3+x2+x1+a1]), (n1-a1)*sizeof(double));
										p+=(n1-a1);
									}							
								}								
							}
	
						}
					}				
				}
			}			
		}
		else
		{
			int a5,a4,a3,a2,a1,i5,i4,i3,i2,i1,n5,n4,n3,n2,n1,x4,x3,x2,x1,y3,y2;
			for(a5=0;a5<r5;a5+=reOrgSize)
			{
				n5 = a5+reOrgSize;
				for(a4=0;a4<r4;a4+=reOrgSize)
				{
					n4 = a4+reOrgSize;
					for(a3=0;a3<r3;a3+=reOrgSize)
					{
						n3 = a3+reOrgSize;
						if(n3 > r3)
							n3 = r3;
						for(a2=0;a2<r2;a2+=reOrgSize)
						{
							n2 = a2+reOrgSize;
							if(n2 > r2)
								n2 = r2;
							for(a1=0;a1<r1;a1+=reOrgSize)
							{
								n1 = a1+reOrgSize;
								if(n1 > r1)
									n1 = r1;
								for(i5=a5;i5<n5;i5++)
								{
									y2 = r2*r1;
									y3 = r3*y2;
									x4 = i5*r4*y3;
									for(i4=a4;i4<n4;i4++)
									{
										x3 = (a4+i4)*y3;
										for(i3=a3;i3<n3;i3++)
										{
											x2 = i3*y2;
											for(i2=a2;i2<n2;i2++)
											{
												x1 = i2*r1;
												memcpy(p, &(oriData[x4+x3+x2+x1+a1]), (n1-a1)*sizeof(double));
												p+=(n1-a1);
											}							
										}								
									}
								}
							}
						}				
					}
				}				
			}
		}
	}

}

void reorganizeData_float_spacefilling(float **newData, float* oriData, int r5, int r4, int r3, int r2, int r1)
{
	if(r2==0)
	{
		*newData = oriData;
	}
	else
	{
		int index = 0;
		int dataLength = computeDataLength(r5,r4,r3,r2,r1);		
		*newData = (float*)malloc(sizeof(float)*dataLength);
		float *p = *newData;
		if(r3==0)
		{
			int a2,a1,i2,i1,n2,n1,x1;
			for(a1=0;a1<r1;a1+=reOrgSize)
			{
				n1 = a1+reOrgSize;
				if(n1 > r1)
					n1 = r1;
				for(a2=0;a2<r2;a2+=reOrgSize)
				{
					n2 = a2+reOrgSize;
					if(n2 > r2)
						n2 = r2;
					for(i2=a2;i2<n2;i2++)
					{
						x1 = i2*r1;
						if(i2%2==0)
						{
							memcpy(p, &(oriData[x1+a1]), (n1-a1)*sizeof(float));
							p+=(n1-a1);
						}
						else
						{
							for(i1=n1-1;i1>=a1;i1--)
								*(p++) = oriData[x1+i1];					
						}
					}			
				}
			}
		}
		else if(r4==0)
		{
			int a3,a2,a1,i3,i2,i1,n3,n2,n1,x2,x1;
			for(a3=0;a3<r3;a3+=reOrgSize)
			{
				n3 = a3+reOrgSize;
				if(n3 > r3)
					n3 = r3;
				for(a1=0;a1<r1;a1+=reOrgSize)
				{
					n1 = a1+reOrgSize;
					if(n1 > r1)
						n1 = r1;
					for(a2=0;a2<r2;a2+=reOrgSize)
					{
						n2 = a2+reOrgSize;
						if(n2 > r2)
							n2 = r2;
						for(i3=a3;i3<n3;i3++)
						{
							x2 = i3*r2*r1;
							for(i2=a2;i2<n2;i2++)
							{
								x1 = i2*r1;
								if(i2%2==0)
								{
									memcpy(p, &(oriData[x2+x1+a1]), (n1-a1)*sizeof(float));
									p+=(n1-a1);
								}
								else
								{
									for(i1=n1-1;i1>=a1;i1--)
										*(p++) = oriData[x2+x1+i1];						
								}
							}							
						}	
					}
				}				
			}	
		}
		else if(r5==0)
		{
			int a4,a3,a2,a1,i4,i3,i2,i1,n4,n3,n2,n1,x3,x2,x1,y2;
			for(a4=0;a4<r4;a4+=reOrgSize)
			{
				n4 = a4+reOrgSize;
				for(a3=0;a3<r3;a3+=reOrgSize)
				{
					n3 = a3+reOrgSize;
					if(n3 > r3)
						n3 = r3;
					for(a1=0;a1<r1;a1+=reOrgSize)
					{
						n1 = a1+reOrgSize;
						if(n1 > r1)
							n1 = r1;
						for(a2=0;a2<r2;a2+=reOrgSize)
						{
							n2 = a2+reOrgSize;
							if(n2 > r2)
								n2 = r2;
							for(i4=a4;i4<n4;i4++)
							{
								y2 = r2*r1;
								x3 = i4*r3*y2;
								for(i3=a3;i3<n3;i3++)
								{
									x2 = i3*y2;
									for(i2=a2;i2<n2;i2++)
									{
										x1 = i2*r1;
										if(i2%2==0)
										{
											memcpy(p, &(oriData[x3+x2+x1+a1]), (n1-a1)*sizeof(float));
											p+=(n1-a1);
										}
										else
										{
											for(i1=n1-1;i1>=a1;i1--)
												*(p++) = oriData[x3+x2+x1+i1];						
										}
									}							
								}								
							}
	
						}
					}				
				}
			}			
		}
		else
		{
			int a5,a4,a3,a2,a1,i5,i4,i3,i2,i1,n5,n4,n3,n2,n1,x4,x3,x2,x1,y3,y2;
			for(a5=0;a5<r5;a5+=reOrgSize)
			{
				n5 = a5+reOrgSize;
				for(a4=0;a4<r4;a4+=reOrgSize)
				{
					n4 = a4+reOrgSize;
					for(a3=0;a3<r3;a3+=reOrgSize)
					{
						n3 = a3+reOrgSize;
						if(n3 > r3)
							n3 = r3;
						for(a1=0;a1<r1;a1+=reOrgSize)
						{
							n1 = a1+reOrgSize;
							if(n1 > r1)
								n1 = r1;
							for(a2=0;a2<r2;a2+=reOrgSize)
							{
								n2 = a2+reOrgSize;
								if(n2 > r2)
									n2 = r2;
								for(i5=a5;i5<n5;i5++)
								{
									y2 = r2*r1;
									y3 = r3*y2;
									x4 = i5*r4*y3;
									for(i4=a4;i4<n4;i4++)
									{
										x3 = (a4+i4)*y3;
										for(i3=a3;i3<n3;i3++)
										{
											x2 = i3*y2;
											for(i2=a2;i2<n2;i2++)
											{
												x1 = i2*r1;
												if(i2%2==0)
												{
													memcpy(p, &(oriData[x4+x3+x2+x1+a1]), (n1-a1)*sizeof(float));
													p+=(n1-a1);
												}
												else
												{
													for(i1=n1-1;i1>=a1;i1--)
														*(p++) = oriData[x4+x3+x2+x1+i1];						
												}
											}								
										}								
									}
								}
							}
						}				
					}
				}				
			}
		}
	}

}

void reorganizeData_double_spacefilling(double **newData, double* oriData, int r5, int r4, int r3, int r2, int r1)
{
	if(r2==0)
	{
		*newData = oriData;
	}
	else
	{
		int index = 0;
		int dataLength = computeDataLength(r5,r4,r3,r2,r1);		
		*newData = (double*)malloc(sizeof(double)*dataLength);
		double *p = *newData;
		if(r3==0)
		{
			int a2,a1,i2,i1,n2,n1,x1;
			for(a1=0;a1<r1;a1+=reOrgSize)
			{
				n1 = a1+reOrgSize;
				if(n1 > r1)
					n1 = r1;
				for(a2=0;a2<r2;a2+=reOrgSize)
				{
					n2 = a2+reOrgSize;
					if(n2 > r2)
						n2 = r2;
					for(i2=a2;i2<n2;i2++)
					{
						x1 = i2*r1;
						if(i2%2==0)
						{
							memcpy(p, &(oriData[x1+a1]), (n1-a1)*sizeof(double));
							p+=(n1-a1);
						}
						else
						{
							for(i1=n1-1;i1>=a1;i1--)
								*(p++) = oriData[x1+i1];					
						}
					}			
				}
			}
		}
		else if(r4==0)
		{
			int a3,a2,a1,i3,i2,i1,n3,n2,n1,x2,x1;
			for(a3=0;a3<r3;a3+=reOrgSize)
			{
				n3 = a3+reOrgSize;
				if(n3 > r3)
					n3 = r3;
				for(a1=0;a1<r1;a1+=reOrgSize)
				{
					n1 = a1+reOrgSize;
					if(n1 > r1)
						n1 = r1;
					for(a2=0;a2<r2;a2+=reOrgSize)
					{
						n2 = a2+reOrgSize;
						if(n2 > r2)
							n2 = r2;
						for(i3=a3;i3<n3;i3++)
						{
							x2 = i3*r2*r1;
							for(i2=a2;i2<n2;i2++)
							{
								x1 = i2*r1;
								if(i2%2==0)
								{
									memcpy(p, &(oriData[x2+x1+a1]), (n1-a1)*sizeof(double));
									p+=(n1-a1);
								}
								else
								{
									for(i1=n1-1;i1>=a1;i1--)
										*(p++) = oriData[x2+x1+i1];						
								}
							}							
						}	
					}
				}				
			}	
		}
		else if(r5==0)
		{
			int a4,a3,a2,a1,i4,i3,i2,i1,n4,n3,n2,n1,x3,x2,x1,y2;
			for(a4=0;a4<r4;a4+=reOrgSize)
			{
				n4 = a4+reOrgSize;
				for(a3=0;a3<r3;a3+=reOrgSize)
				{
					n3 = a3+reOrgSize;
					if(n3 > r3)
						n3 = r3;
					for(a1=0;a1<r1;a1+=reOrgSize)
					{
						n1 = a1+reOrgSize;
						if(n1 > r1)
							n1 = r1;
						for(a2=0;a2<r2;a2+=reOrgSize)
						{
							n2 = a2+reOrgSize;
							if(n2 > r2)
								n2 = r2;
							for(i4=a4;i4<n4;i4++)
							{
								y2 = r2*r1;
								x3 = i4*r3*y2;
								for(i3=a3;i3<n3;i3++)
								{
									x2 = i3*y2;
									for(i2=a2;i2<n2;i2++)
									{
										x1 = i2*r1;
										if(i2%2==0)
										{
											memcpy(p, &(oriData[x3+x2+x1+a1]), (n1-a1)*sizeof(double));
											p+=(n1-a1);
										}
										else
										{
											for(i1=n1-1;i1>=a1;i1--)
												*(p++) = oriData[x3+x2+x1+i1];						
										}
									}							
								}								
							}
	
						}
					}				
				}
			}			
		}
		else
		{
			int a5,a4,a3,a2,a1,i5,i4,i3,i2,i1,n5,n4,n3,n2,n1,x4,x3,x2,x1,y3,y2;
			for(a5=0;a5<r5;a5+=reOrgSize)
			{
				n5 = a5+reOrgSize;
				for(a4=0;a4<r4;a4+=reOrgSize)
				{
					n4 = a4+reOrgSize;
					for(a3=0;a3<r3;a3+=reOrgSize)
					{
						n3 = a3+reOrgSize;
						if(n3 > r3)
							n3 = r3;
						for(a1=0;a1<r1;a1+=reOrgSize)
						{
							n1 = a1+reOrgSize;
							if(n1 > r1)
								n1 = r1;
							for(a2=0;a2<r2;a2+=reOrgSize)
							{
								n2 = a2+reOrgSize;
								if(n2 > r2)
									n2 = r2;
								for(i5=a5;i5<n5;i5++)
								{
									y2 = r2*r1;
									y3 = r3*y2;
									x4 = i5*r4*y3;
									for(i4=a4;i4<n4;i4++)
									{
										x3 = (a4+i4)*y3;
										for(i3=a3;i3<n3;i3++)
										{
											x2 = i3*y2;
											for(i2=a2;i2<n2;i2++)
											{
												x1 = i2*r1;
												if(i2%2==0)
												{
													memcpy(p, &(oriData[x4+x3+x2+x1+a1]), (n1-a1)*sizeof(double));
													p+=(n1-a1);
												}
												else
												{
													for(i1=n1-1;i1>=a1;i1--)
														*(p++) = oriData[x4+x3+x2+x1+i1];						
												}
											}								
										}								
									}
								}
							}
						}				
					}
				}				
			}
		}
	}

}
