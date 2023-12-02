/* these are extra and sort of hard to produce automatically */
void assemble_ap_matrix_1(A_MATRIX_PREC *,int ,int ,my_boolean *,my_boolean *,int ,int ,
			  int *,int *,struct flt *,struct med *);
void assemble_ap_matrix_2(A_MATRIX_PREC *,int ,int ,my_boolean *,my_boolean *,int ,int ,
			  int *,int *,struct flt *,struct med *);
void assemble_ap_matrix_3(A_MATRIX_PREC *,int ,int ,my_boolean *,my_boolean *,int ,int ,
			  int *,int *,struct flt *,struct med *);
void assemble_ap_matrix_4(A_MATRIX_PREC *,int ,int ,my_boolean *,my_boolean *,int ,int ,
			  int *,int *,struct flt *,struct med *);

int par_assemble_a_matrix(int ,my_boolean *,int ,int *,struct flt *,struct med *);

void add_quake_stress_1(my_boolean *,COMP_PRECISION *,int ,struct flt *,struct med *);
void add_quake_stress_2(my_boolean *,COMP_PRECISION *,int ,struct flt *,struct med *);
void add_quake_stress_3(my_boolean *,COMP_PRECISION *,int ,struct flt *,struct med *);
void add_quake_stress_4(my_boolean *,COMP_PRECISION *,int ,struct flt *,struct med *);
my_boolean check_coulomb_stress_feedback_1(int ,int,struct flt *,struct med *,my_boolean ,my_boolean,int *,COMP_PRECISION);
my_boolean check_coulomb_stress_feedback_2(int ,int,struct flt *,struct med *,my_boolean ,my_boolean,int *,COMP_PRECISION);
my_boolean check_coulomb_stress_feedback_3(int ,int,struct flt *,struct med *,my_boolean ,my_boolean,int *,COMP_PRECISION);
my_boolean check_coulomb_stress_feedback_4(int ,int,struct flt *,struct med *,my_boolean ,my_boolean,int *,COMP_PRECISION);
void assemble_a_matrix_1(A_MATRIX_PREC *,int ,my_boolean *,int ,int *,struct flt *,struct med *);
void assemble_a_matrix_2(A_MATRIX_PREC *,int ,my_boolean *,int ,int *,struct flt *,struct med *);
void assemble_a_matrix_3(A_MATRIX_PREC *,int ,my_boolean *,int ,int *,struct flt *,struct med *);
void assemble_a_matrix_4(A_MATRIX_PREC *,int ,my_boolean *,int ,int *,struct flt *,struct med *);
