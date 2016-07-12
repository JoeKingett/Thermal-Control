function [n1,n2,n3,n4,n5,n6]=normal(acend_node,i,w,true_anomoly)

%convert to radians
acend_node=acend_node*pi/180;
i=i*pi/180;
w=w*pi/180;

norm1_local=[1; 0; 0];
norm2_local=[-1; 0; 0];
norm3_local=[0; 1; 0];
norm4_local=[0; -1; 0];
norm5_local=[0; 0; 1];
norm6_local=[0; 0; -1];

T=[cos(acend_node)*cos(w+true_anomoly)-sin(acend_node)*cos(i)*sin(w+true_anomoly)...
              -cos(acend_node)*sin(w+true_anomoly)-sin(acend_node)*cos(i)*cos(w+true_anomoly)...
              sin(acend_node)*sin(i);...
   sin(acend_node)*cos(w+true_anomoly)+cos(acend_node)*cos(i)*sin(w+true_anomoly)...
              -sin(acend_node)*sin(w+true_anomoly)+cos(acend_node)*cos(i)*cos(w+true_anomoly)...
              -cos(acend_node)*sin(i);...
  sin(i)*sin(w+true_anomoly)     sin(i)*cos(w+true_anomoly)        cos(i)];

n1=T*norm1_local;
n2=T*norm2_local;
n3=T*norm3_local;
n4=T*norm4_local;
n5=T*norm5_local;
n6=T*norm6_local;


   