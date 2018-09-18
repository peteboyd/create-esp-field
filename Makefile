FC = gfortran
FFLAGS = -O3 
all: create-esp-field 

create-esp-field: 
	$(FC) $(FFLAGS) create_esp_field.f -o create-esp-field 

clean:
	rm -f create-esp-field 
