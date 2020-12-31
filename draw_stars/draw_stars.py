# Kendra Ludwig (kel334@nau.edu)


import turtle


def main():
    # create the turtle and settings
    keya = turtle.Turtle()
    turtle.bgcolor('black')
    turtle.screensize(500, 500)
    turtle.speed(0)
    turtle.tracer(1, 0)
    turtle.forward(10)
    turtle.color('white')

    # open the input file and get x, y, and m values
    with open('stars.txt', 'r') as star_file:
        for line in star_file:
            value = line.split()
            x = float(value[0])
            y = float(value[1])
            m = float(value[4])
            side = round(10/(m + 2))

            # restrict the size of the stars
            if side > 10:
                side = 10
            if side < 1:
                side = 1

            # draw the stars
            turtle.penup()
            turtle.goto(x*250, y*250)
            turtle.pendown()
            turtle.begin_fill()
            draw_square(side)
            turtle.end_fill()

    # call the functions for the constellations
    turtle.pencolor('yellow')
    big_dipper()
    ursa_major()
    ursa_minor()
    gemini()
    bootes()
    cassiopeia()
    cygnet()
    hydra()

    turtle.ht()
    turtle.exitonclick()


# draw the square for reference
def draw_square(a):
    for i in range(4):
        turtle.forward(a)
        turtle.right(90)


# the constellations begin here: the big dipper
def big_dipper():
    with open('stars.txt', 'r') as original_star_file:
        the_big_dipper = original_star_file.readlines()
    big_dipper_stars = []

    # place all the necessary lines into a new file
    for stars in the_big_dipper:
        if 'BENETNASCH' in stars:
            big_dipper_stars.append(stars)
        if 'MIZAR' in stars:
            big_dipper_stars.append(stars)
        if 'ALIOTH' in stars:
            big_dipper_stars.append(stars)
        if 'MEGREZ' in stars:
            big_dipper_stars.append(stars)
        if 'PHECDA' in stars:
            big_dipper_stars.append(stars)
        if 'MERAK' in stars:
            big_dipper_stars.append(stars)
        if 'DUBHE' in stars:
            big_dipper_stars.append(stars)

    with open('big_dipper_stars.txt', 'w') as bigdipper_stars:
        bigdipper_stars.writelines(big_dipper_stars)

    # using the new file, assign values to the coordinates
    with open('big_dipper_stars.txt', 'r') as stars_big_dipper:
        for lines in stars_big_dipper:
            coordinates = lines.split()
            if 'BENETNASCH' in lines:
                benetnasch_x = coordinates[0]
                benetnasch_y = coordinates[1]
            if 'MIZAR' in lines:
                mizar_x = coordinates[0]
                mizar_y = coordinates[1]
            if 'ALIOTH' in lines:
                alioth_x = coordinates[0]
                alioth_y = coordinates[1]
            if 'MEGREZ' in lines:
                megrez_x = coordinates[0]
                megrez_y = coordinates[1]
            if 'PHECDA' in lines:
                phecda_x = coordinates[0]
                phecda_y = coordinates[1]
            if 'MERAK' in lines:
                merak_x = coordinates[0]
                merak_y = coordinates[1]
            if 'DUBHE' in lines:
                dubhe_x = coordinates[0]
                dubhe_y = coordinates[1]

        # draw the lines between the stars
        turtle.penup()
        turtle.goto(float(benetnasch_x)*250, float(benetnasch_y)*250)
        turtle.pendown()
        turtle.goto(float(mizar_x)*250, float(mizar_y)*250)
        turtle.goto(float(alioth_x)*250, float(alioth_y)*250)
        turtle.goto(float(megrez_x)*250, float(megrez_y)*250)
        turtle.goto(float(phecda_x)*250, float(phecda_y)*250)
        turtle.goto(float(merak_x)*250, float(merak_y)*250)
        turtle.goto(float(dubhe_x)*250, float(dubhe_y)*250)
        turtle.goto(float(megrez_x)*250, float(megrez_y)*250)


# repeat process above from big dipper with other constellations
def ursa_major():
    with open('stars.txt', 'r') as original_star_file:
        ursa_major = original_star_file.readlines()

    ursa_major_stars = []
    for stars in ursa_major:
        if 'BENETNASCH' in stars:
            ursa_major_stars.append(stars)
        if 'MIZAR' in stars:
            ursa_major_stars.append(stars)
        if 'ALIOTH' in stars:
            ursa_major_stars.append(stars)
        if 'MEGREZ' in stars:
            ursa_major_stars.append(stars)
        if 'PHECDA' in stars:
            ursa_major_stars.append(stars)
        if 'MERAK' in stars:
            ursa_major_stars.append(stars)
        if 'DUBHE' in stars:
            ursa_major_stars.append(stars)
        if 'UMA OMICRON' in stars:
            ursa_major_stars.append(stars)
        if 'UMA PHI' in stars:
            ursa_major_stars.append(stars)
        if 'UMA UPSILON' in stars:
            ursa_major_stars.append(stars)
        if 'UMA LAMBDA' in stars:
            ursa_major_stars.append(stars)
        if 'TALITHA' in stars:
            ursa_major_stars.append(stars)
        if 'UMA KAPPA' in stars:
            ursa_major_stars.append(stars)
        if 'UMA MU' in stars:
            ursa_major_stars.append(stars)
        if 'UMA NU' in stars:
            ursa_major_stars.append(stars)
        if 'UMA PSI' in stars:
            ursa_major_stars.append(stars)
        if 'UMA CHI' in stars:
            ursa_major_stars.append(stars)

    with open('ursa_major_stars.txt', 'w') as ursamajor_stars:
        ursamajor_stars.writelines(ursa_major_stars)

    with open('ursa_major_stars.txt', 'r') as stars_ursa_major:
        for lines in stars_ursa_major:
            coordinates = lines.split()
            if 'BENETNASCH' in lines:
                benetnasch_x = coordinates[0]
                benetnasch_y = coordinates[1]
            if 'MIZAR' in lines:
                mizar_x = coordinates[0]
                mizar_y = coordinates[1]
            if 'ALIOTH' in lines:
                alioth_x = coordinates[0]
                alioth_y = coordinates[1]
            if 'MEGREZ' in lines:
                megrez_x = coordinates[0]
                megrez_y = coordinates[1]
            if 'PHECDA' in lines:
                phecda_x = coordinates[0]
                phecda_y = coordinates[1]
            if 'MERAK' in lines:
                merak_x = coordinates[0]
                merak_y = coordinates[1]
            if 'DUBHE' in lines:
                dubhe_x = coordinates[0]
                dubhe_y = coordinates[1]
            if 'UMA OMICRON' in lines:
                uma_omicron_x = coordinates[0]
                uma_omicron_y = coordinates[1]
            if 'UMA PHI' in lines:
                uma_phi_x = coordinates[0]
                uma_phi_y = coordinates[1]
            if 'UMA UPSILON' in lines:
                uma_upsilon_x = coordinates[0]
                uma_upsilon_y = coordinates[1]
            if 'UMA LAMBDA' in lines:
                uma_lambda_x = coordinates[0]
                uma_lambda_y = coordinates[1]
            if 'TALITHA' in lines:
                talitha_x = coordinates[0]
                talitha_y = coordinates[1]
            if 'UMA KAPPA' in lines:
                uma_kappa_x = coordinates[0]
                uma_kappa_y = coordinates[1]
            if 'UMA MU' in lines:
                uma_mu_x = coordinates[0]
                uma_mu_y = coordinates[1]
            if 'UMA NU' in lines:
                uma_nu_x = coordinates[0]
                uma_nu_y = coordinates[1]
            if 'UMA PSI' in lines:
                uma_psi_x = coordinates[0]
                uma_psi_y = coordinates[1]
            if 'UMA CHI' in lines:
                uma_chi_x = coordinates[0]
                uma_chi_y = coordinates[1]

        turtle.penup()
        turtle.goto(float(benetnasch_x)*250, float(benetnasch_y)*250)
        turtle.pendown()
        turtle.goto(float(mizar_x)*250, float(mizar_y)*250)
        turtle.goto(float(alioth_x)*250, float(alioth_y)*250)
        turtle.goto(float(megrez_x)*250, float(megrez_y)*250)
        turtle.goto(float(phecda_x)*250, float(phecda_y)*250)
        turtle.goto(float(merak_x)*250, float(merak_y)*250)
        turtle.goto(float(dubhe_x)*250, float(dubhe_y)*250)
        turtle.goto(float(megrez_x)*250, float(megrez_y)*250)
        turtle.goto(float(dubhe_x)*250, float(dubhe_y)*250)
        turtle.goto(float(uma_omicron_x)*250, float(uma_omicron_y)*250)
        turtle.goto(float(talitha_x)*250, float(talitha_y)*250)
        turtle.penup()
        turtle.goto(float(uma_lambda_x)*250, float(uma_lambda_y)*250)
        turtle.pendown()
        turtle.goto(float(uma_phi_x)*250, float(uma_phi_y)*250)
        turtle.goto(float(uma_upsilon_x)*250, float(uma_upsilon_y)*250)
        turtle.penup()
        turtle.goto(float(uma_kappa_x)*250, float(uma_kappa_y)*250)
        turtle.pendown()
        turtle.goto(float(uma_mu_x)*250, float(uma_mu_y)*250)
        turtle.goto(float(uma_nu_x)*250, float(uma_nu_y)*250)
        turtle.goto(float(uma_psi_x)*250, float(uma_psi_y)*250)
        turtle.goto(float(uma_chi_x)*250, float(uma_chi_y)*250)
        turtle.goto(float(benetnasch_x)*250, float(benetnasch_y)*250)


def ursa_minor():
    with open('stars.txt', 'r') as original_star_file:
        ursa_minor = original_star_file.readlines()

    ursa_minor_stars = []
    for stars in ursa_minor:
        if 'POLARIS' in stars:
            ursa_minor_stars.append(stars)
        if 'UMI DELTA' in stars:
            ursa_minor_stars.append(stars)
        if 'UMI EPSILON' in stars:
            ursa_minor_stars.append(stars)
        if 'UMI ZETA' in stars:
            ursa_minor_stars.append(stars)
        if 'KOCHAB' in stars:
            ursa_minor_stars.append(stars)
        if 'UMI GAMMA' in stars:
            ursa_minor_stars.append(stars)
        if 'UMI ETA' in stars:
            ursa_minor_stars.append(stars)

    with open('ursa_minor_stars.txt', 'w') as ursa_minor_star_file:
        ursa_minor_star_file.writelines(ursa_minor_stars)

    with open('ursa_minor_stars.txt', 'r') as stars_ursa_minor:
        for lines in stars_ursa_minor:
            coordinates = lines.split()
            if 'POLARIS' in lines:
                polaris_x = coordinates[0]
                polaris_y = coordinates[1]
            if 'UMI DELTA' in lines:
                umi_delta_x = coordinates[0]
                umi_delta_y = coordinates[1]
            if 'UMI EPSILON' in lines:
                umi_epsilon_x = coordinates[0]
                umi_epsilon_y = coordinates[1]
            if 'UMI ZETA' in lines:
                umi_zeta_x = coordinates[0]
                umi_zeta_y = coordinates[1]
            if 'KOCHAB' in lines:
                kochab_x = coordinates[0]
                kochab_y = coordinates[1]
            if 'UMI GAMMA' in lines:
                umi_gamma_x = coordinates[0]
                umi_gamma_y = coordinates[1]
            if 'UMI ETA' in lines:
                umi_eta_x = coordinates[0]
                umi_eta_y = coordinates[1]

        turtle.penup()
        turtle.goto(float(polaris_x)*250, float(polaris_y)*250)
        turtle.pendown()
        turtle.goto(float(umi_delta_x)*250, float(umi_delta_y)*250)
        turtle.goto(float(umi_epsilon_x)*250, float(umi_epsilon_y)*250)
        turtle.goto(float(umi_zeta_x)*250, float(umi_zeta_y)*250)
        turtle.goto(float(kochab_x)*250, float(kochab_y)*250)
        turtle.goto(float(umi_gamma_x)*250, float(umi_gamma_y)*250)
        turtle.goto(float(umi_eta_x)*250, float(umi_eta_y)*250)
        turtle.goto(float(umi_zeta_x)*250, float(umi_zeta_y)*250)


def gemini():
    with open('stars.txt', 'r') as original_star_file:
        gemini = original_star_file.readlines()

    gemini_stars = []
    for stars in gemini:
        if 'POLLUX' in stars:
            gemini_stars.append(stars)
        if 'GEM UPSILON' in stars:
            gemini_stars.append(stars)
        if 'GEM KAPPA' in stars:
            gemini_stars.append(stars)
        if 'GEM IOTA' in stars:
            gemini_stars.append(stars)
        if 'WASAT' in stars:
            gemini_stars.append(stars)
        if 'GEM LAMBDA' in stars:
            gemini_stars.append(stars)
        if 'GEM ZETA' in stars:
            gemini_stars.append(stars)
        if 'GEM GAMMA' in stars:
            gemini_stars.append(stars)
        if 'GEM XI' in stars:
            gemini_stars.append(stars)
        if 'GEM RHO' in stars:
            gemini_stars.append(stars)
        if 'CASTOR' in stars:
            gemini_stars.append(stars)
        if 'GEM TAU' in stars:
            gemini_stars.append(stars)
        if 'GEM THETA' in stars:
            gemini_stars.append(stars)
        if 'MEBSUTA' in stars:
            gemini_stars.append(stars)
        if 'GEM MU' in stars:
            gemini_stars.append(stars)
        if 'GEM ETA' in stars:
            gemini_stars.append(stars)

    with open('gemini_stars.txt', 'w') as gemini_star_file:
        gemini_star_file.writelines(gemini_stars)

    with open('gemini_stars.txt', 'r') as stars_gemini:
        for lines in stars_gemini:
            coordinates = lines.split()
            if 'POLLUX' in lines:
                pollux_x = coordinates[0]
                pollux_y = coordinates[1]
            if 'GEM UPSILON' in lines:
                gem_upsilon_x = coordinates[0]
                gem_upsilon_y = coordinates[1]
            if 'GEM KAPPA' in lines:
                gem_kappa_x = coordinates[0]
                gem_kappa_y = coordinates[1]
            if 'GEM IOTA' in lines:
                gem_iota_x = coordinates[0]
                gem_iota_y = coordinates[1]
            if 'WASAT' in lines:
                wasat_x = coordinates[0]
                wasat_y = coordinates[1]
            if 'GEM LAMBDA' in lines:
                gem_lambda_x = coordinates[0]
                gem_lambda_y = coordinates[1]
            if 'GEM ZETA' in lines:
                gem_zeta_x = coordinates[0]
                gem_zeta_y = coordinates[1]
            if 'GEM GAMMA' in lines:
                gem_gamma_x = coordinates[0]
                gem_gamma_y = coordinates[1]
            if 'GEM XI' in lines:
                gem_xi_x = coordinates[0]
                gem_xi_y = coordinates[1]
            if 'GEM RHO' in lines:
                gem_rho_x = coordinates[0]
                gem_rho_y = coordinates[1]
            if 'CASTOR' in lines:
                castor_x = coordinates[0]
                castor_y = coordinates[1]
            if 'GEM TAU' in lines:
                gem_tau_x = coordinates[0]
                gem_tau_y = coordinates[1]
            if 'GEM THETA' in lines:
                gem_theta_x = coordinates[0]
                gem_theta_y = coordinates[1]
            if 'MEBSUTA' in lines:
                mebsuta_x = coordinates[0]
                mebsuta_y = coordinates[1]
            if 'GEM MU' in lines:
                gem_mu_x = coordinates[0]
                gem_mu_y = coordinates[1]
            if 'GEM ETA' in lines:
                gem_eta_x = coordinates[0]
                gem_eta_y = coordinates[1]

        turtle.penup()
        turtle.goto(float(pollux_x)*250, float(pollux_y)*250)
        turtle.pendown()
        turtle.goto(float(gem_upsilon_x)*250, float(gem_upsilon_y)*250)
        turtle.goto(float(gem_kappa_x)*250, float(gem_kappa_y)*250)
        turtle.goto(float(gem_upsilon_x)*250, float(gem_upsilon_y)*250)
        turtle.goto(float(gem_iota_x)*250, float(gem_iota_y)*250)
        turtle.goto(float(gem_upsilon_x)*250, float(gem_upsilon_y)*250)
        turtle.goto(float(wasat_x)*250, float(wasat_y)*250)
        turtle.goto(float(gem_lambda_x)*250, float(gem_lambda_y)*250)
        turtle.goto(float(gem_xi_x)*250, float(gem_xi_y)*250)
        turtle.penup()
        turtle.goto(float(wasat_x)*250, float(wasat_y)*250)
        turtle.pendown()
        turtle.goto(float(gem_zeta_x)*250, float(gem_zeta_y)*250)
        turtle.goto(float(gem_gamma_x)*250, float(gem_gamma_y)*250)
        turtle.penup()
        turtle.goto(float(castor_x)*250, float(castor_y)*250)
        turtle.pendown()
        turtle.goto(float(gem_rho_x)*250, float(gem_rho_y)*250)
        turtle.goto(float(gem_tau_x)*250, float(gem_tau_y)*250)
        turtle.goto(float(gem_iota_x)*250, float(gem_iota_y)*250)
        turtle.goto(float(gem_tau_x)*250, float(gem_tau_y)*250)
        turtle.goto(float(gem_theta_x)*250, float(gem_theta_y)*250)
        turtle.goto(float(gem_tau_x)*250, float(gem_tau_y)*250)
        turtle.goto(float(mebsuta_x)*250, float(mebsuta_y)*250)
        turtle.goto(float(gem_mu_x)*250, float(gem_mu_y)*250)
        turtle.goto(float(gem_eta_x)*250, float(gem_eta_y)*250)


def bootes():
    with open('stars.txt', 'r') as original_star_file:
        bootes = original_star_file.readlines()

    bootes_stars = []
    for stars in bootes:
        if 'ARCTURUS' in stars:
            bootes_stars.append(stars)
        if 'BOO PI' in stars:
            bootes_stars.append(stars)
        if 'PULCHERRIMA' in stars:
            bootes_stars.append(stars)
        if 'MUPHRID' in stars:
            bootes_stars.append(stars)
        if 'BOO UPSILON' in stars:
            bootes_stars.append(stars)
        if 'BOO TAU' in stars:
            bootes_stars.append(stars)
        if 'NEKKAR' in stars:
            bootes_stars.append(stars)
        if 'BOO DELTA' in stars:
            bootes_stars.append(stars)
        if 'BOO MU' in stars:
            bootes_stars.append(stars)
        if 'BOO GAMMA' in stars:
            bootes_stars.append(stars)
        if 'BOO SIGMA' in stars:
            bootes_stars.append(stars)
        if 'BOO RHO' in stars:
            bootes_stars.append(stars)
        if 'BOO THETA' in stars:
            bootes_stars.append(stars)
        if 'BOO KAPPA' in stars:
            bootes_stars.append(stars)
        if 'BOO LAMBDA' in stars:
            bootes_stars.append(stars)

    with open('bootes_stars.txt', 'w') as bootes_star_file:
        bootes_star_file.writelines(bootes_stars)

    with open('bootes_stars.txt', 'r') as stars_bootes:
        for lines in stars_bootes:
            coordinates = lines.split()
            if 'ARCTURUS' in lines:
                arcturus_x = coordinates[0]
                arcturus_y = coordinates[1]
            if 'BOO PI' in lines:
                boo_pi_x = coordinates[0]
                boo_pi_y = coordinates[1]
            if 'PULCHERRIMA' in lines:
                pulcherrima_x = coordinates[0]
                pulcherrima_y = coordinates[1]
            if 'MUPHRID' in lines:
                muphrid_x = coordinates[0]
                muphrid_y = coordinates[1]
            if 'BOO UPSILON' in lines:
                boo_upsilon_x = coordinates[0]
                boo_upsilon_y = coordinates[1]
            if 'BOO TAU' in lines:
                boo_tau_x = coordinates[0]
                boo_tau_y = coordinates[1]
            if 'NEKKAR' in lines:
                nekkar_x = coordinates[0]
                nekkar_y = coordinates[1]
            if 'BOO DELTA' in lines:
                boo_delta_x = coordinates[0]
                boo_delta_y = coordinates[1]
            if 'BOO MU' in lines:
                boo_mu_x = coordinates[0]
                boo_mu_y = coordinates[1]
            if 'BOO GAMMA' in lines:
                boo_gamma_x = coordinates[0]
                boo_gamma_y = coordinates[1]
            if 'BOO SIGMA' in lines:
                boo_sigma_x = coordinates[0]
                boo_sigma_y = coordinates[1]
            if 'BOO RHO' in lines:
                boo_rho_x = coordinates[0]
                boo_rho_y = coordinates[1]
            if 'BOO THETA' in lines:
                boo_theta_x = coordinates[0]
                boo_theta_y = coordinates[1]
            if 'BOO KAPPA' in lines:
                boo_kappa_x = coordinates[0]
                boo_kappa_y = coordinates[1]
            if 'BOO LAMBDA' in lines:
                boo_lambda_x = coordinates[0]
                boo_lambda_y = coordinates[1]

        turtle.penup()
        turtle.goto(float(arcturus_x)*250, float(arcturus_y)*250)
        turtle.pendown()
        turtle.goto(float(pulcherrima_x)*250, float(pulcherrima_y)*250)
        turtle.goto(float(boo_delta_x)*250, float(boo_delta_y)*250)
        turtle.goto(float(boo_mu_x)*250, float(boo_mu_y)*250)
        turtle.goto(float(nekkar_x)*250, float(nekkar_y)*250)
        turtle.goto(float(boo_gamma_x)*250, float(boo_gamma_y)*250)
        turtle.goto(float(boo_lambda_x)*250, float(boo_lambda_y)*250)
        turtle.goto(float(boo_kappa_x)*250, float(boo_kappa_y)*250)
        turtle.goto(float(boo_theta_x)*250, float(boo_theta_y)*250)
        turtle.goto(float(boo_lambda_x)*250, float(boo_lambda_y)*250)
        turtle.goto(float(boo_gamma_x)*250, float(boo_gamma_y)*250)
        turtle.goto(float(boo_rho_x)*250, float(boo_rho_y)*250)
        turtle.goto(float(boo_sigma_x)*250, float(boo_sigma_y)*250)
        turtle.goto(float(pulcherrima_x)*250, float(pulcherrima_y)*250)
        turtle.goto(float(arcturus_x)*250, float(arcturus_y)*250)
        turtle.goto(float(muphrid_x)*250, float(muphrid_y)*250)
        turtle.goto(float(boo_upsilon_x)*250, float(boo_upsilon_y)*250)
        turtle.goto(float(boo_tau_x)*250, float(boo_tau_y)*250)
        turtle.penup()
        turtle.goto(float(arcturus_x)*250, float(arcturus_y)*250)
        turtle.pendown()
        turtle.goto(float(boo_pi_x)*250, float(boo_pi_y)*250)
        turtle.goto(float(pulcherrima_x)*250, float(pulcherrima_y)*250)
        turtle.penup()
        turtle.goto(float(nekkar_x)*250, float(nekkar_y)*250)
        turtle.pendown()
        turtle.goto(float(boo_delta_x)*250, float(boo_delta_y)*250)


def cassiopeia():
    with open('stars.txt', 'r') as original_star_file:
        cassiopeia = original_star_file.readlines()
    cassiopeia_stars = []

    for stars in cassiopeia:
        if 'CAS EPSILON' in stars:
            cassiopeia_stars.append(stars)
        if 'CAS DELTA' in stars:
            cassiopeia_stars.append(stars)
        if 'CAS GAMMA' in stars:
            cassiopeia_stars.append(stars)
        if 'SCHEDAR' in stars:
            cassiopeia_stars.append(stars)
        if 'CAPH' in stars:
            cassiopeia_stars.append(stars)
        if 'CAS KAPPA' in stars:
            cassiopeia_stars.append(stars)
        if 'CAS IOTA' in stars:
            cassiopeia_stars.append(stars)
        if 'CAS ZETA' in stars:
            cassiopeia_stars.append(stars)
        if 'CAS THETA' in stars:
            cassiopeia_stars.append(stars)

    with open('cassiopeia_stars.txt', 'w') as cassio_stars:
        cassio_stars.writelines(cassiopeia_stars)

    with open('cassiopeia_stars.txt', 'r') as stars_cassiopeia:
        for lines in stars_cassiopeia:
            coordinates = lines.split()
            if 'CAS EPSILON' in lines:
                cas_epsilon_x = coordinates[0]
                cas_epsilon_y = coordinates[1]
            if 'CAS DELTA' in lines:
                cas_delta_x = coordinates[0]
                cas_delta_y = coordinates[1]
            if 'CAS GAMMA' in lines:
                cas_gamma_x = coordinates[0]
                cas_gamma_y = coordinates[1]
            if 'SCHEDAR' in lines:
                schedar_x = coordinates[0]
                schedar_y = coordinates[1]
            if 'CAPH' in lines:
                caph_x = coordinates[0]
                caph_y = coordinates[1]
            if 'CAS KAPPA' in lines:
                cas_kappa_x = coordinates[0]
                cas_kappa_y = coordinates[1]
            if 'CAS IOTA' in lines:
                cas_iota_x = coordinates[0]
                cas_iota_y = coordinates[1]
            if 'CAS ZETA' in lines:
                cas_zeta_x = coordinates[0]
                cas_zeta_y = coordinates[1]
            if 'CAS THETA' in lines:
                cas_theta_x = coordinates[0]
                cas_theta_y = coordinates[1]

        turtle.penup()
        turtle.goto(float(cas_epsilon_x)*250, float(cas_epsilon_y)*250)
        turtle.pendown()
        turtle.goto(float(cas_delta_x)*250, float(cas_delta_y)*250)
        turtle.goto(float(cas_gamma_x)*250, float(cas_gamma_y)*250)
        turtle.goto(float(schedar_x)*250, float(schedar_y)*250)
        turtle.goto(float(caph_x)*250, float(caph_y)*250)
        turtle.penup()
        turtle.goto(float(cas_gamma_x)*250, float(cas_gamma_y)*250)
        turtle.pendown()
        turtle.goto(float(cas_kappa_x)*250, float(cas_kappa_y)*250)
        turtle.goto(float(cas_iota_x)*250, float(cas_iota_y)*250)
        turtle.penup()
        turtle.goto(float(schedar_x)*250, float(schedar_y)*250)
        turtle.pendown()
        turtle.goto(float(cas_zeta_x)*250, float(cas_zeta_y)*250)
        turtle.goto(float(cas_theta_x)*250, float(cas_theta_y)*250)


def cygnet():
    with open('stars.txt', 'r') as original_star_file:
        cygnet = original_star_file.readlines()
    cygnet_stars = []

    for stars in cygnet:
        if 'CYG CHI' in stars:
            cygnet_stars.append(stars)
        if 'ALBIREO' in stars:
            cygnet_stars.append(stars)
        if 'CYG ETA' in stars:
            cygnet_stars.append(stars)
        if 'CYG GAMMA' in stars:
            cygnet_stars.append(stars)
        # renamed deneb to cyg alpha to get right lines
        if 'CYG ALPHA' in stars:
            cygnet_stars.append(stars)
        if 'CYG DELTA' in stars:
            cygnet_stars.append(stars)
        if 'CYG EPSILON' in stars:
            cygnet_stars.append(stars)
        if 'CYG KAPPA' in stars:
            cygnet_stars.append(stars)
        if 'CYG IOTA' in stars:
            cygnet_stars.append(stars)
        if 'CYG ZETA' in stars:
            cygnet_stars.append(stars)
        if 'CYG NU' in stars:
            cygnet_stars.append(stars)
        if 'CYG XI' in stars:
            cygnet_stars.append(stars)
        if 'CYG SIGMA' in stars:
            cygnet_stars.append(stars)
        if 'CYG TAU' in stars:
            cygnet_stars.append(stars)

    with open('cygnet_stars.txt', 'w') as cygn_stars:
        cygn_stars.writelines(cygnet_stars)

    with open('cygnet_stars.txt', 'r') as stars_cygnet:
        for lines in stars_cygnet:
            coordinates = lines.split()
            if 'CYG CHI' in lines:
                cyg_chi_x = coordinates[0]
                cyg_chi_y = coordinates[1]
            if 'ALBIREO' in lines:
                albireo_x = coordinates[0]
                albireo_y = coordinates[1]
            if 'CYG ETA' in lines:
                cyg_eta_x = coordinates[0]
                cyg_eta_y = coordinates[1]
            if 'CYG GAMMA' in lines:
                cyg_gamma_x = coordinates[0]
                cyg_gamma_y = coordinates[1]
            if 'DENEB' in lines:
                deneb_x = coordinates[0]
                deneb_y = coordinates[1]
            if 'CYG DELTA' in lines:
                cyg_delta_x = coordinates[0]
                cyg_delta_y = coordinates[1]
            if 'CYG EPSILON' in lines:
                cyg_epsilon_x = coordinates[0]
                cyg_epsilon_y = coordinates[1]
            if 'CYG KAPPA' in lines:
                cyg_kappa_x = coordinates[0]
                cyg_kappa_y = coordinates[1]
            if 'CYG IOTA' in lines:
                cyg_iota_x = coordinates[0]
                cyg_iota_y = coordinates[1]
            if 'CYG ZETA' in lines:
                cyg_zeta_x = coordinates[0]
                cyg_zeta_y = coordinates[1]
            if 'CYG NU' in lines:
                cyg_nu_x = coordinates[0]
                cyg_nu_y = coordinates[1]
            if 'CYG XI' in lines:
                cyg_xi_x = coordinates[0]
                cyg_xi_y = coordinates[1]
            if 'CYG SIGMA' in lines:
                cyg_sigma_x = coordinates[0]
                cyg_sigma_y = coordinates[1]
            if 'CYG TAU' in lines:
                cyg_tau_x = coordinates[0]
                cyg_tau_y = coordinates[1]

        turtle.penup()
        turtle.goto(float(albireo_x)*250, float(albireo_y)*250)
        turtle.pendown()
        turtle.goto(float(cyg_chi_x)*250, float(cyg_chi_y)*250)
        turtle.goto(float(cyg_eta_x)*250, float(cyg_eta_y)*250)
        turtle.goto(float(cyg_gamma_x)*250, float(cyg_gamma_y)*250)
        turtle.goto(float(cyg_delta_x)*250, float(cyg_delta_y)*250)
        turtle.goto(float(cyg_kappa_x)*250, float(cyg_kappa_y)*250)
        turtle.goto(float(cyg_iota_x)*250, float(cyg_iota_y)*250)
        turtle.goto(float(deneb_x)*250, float(deneb_y)*250)
        turtle.goto(float(cyg_nu_x)*250, float(cyg_nu_y)*250)
        turtle.goto(float(cyg_sigma_x)*250, float(cyg_sigma_y)*250)
        turtle.goto(float(cyg_tau_x)*250, float(cyg_tau_y)*250)
        turtle.penup()
        turtle.goto(float(deneb_x)*250, float(deneb_y)*250)
        turtle.pendown()
        turtle.goto(float(cyg_gamma_x)*250, float(cyg_gamma_y)*250)
        turtle.goto(float(cyg_epsilon_x)*250, float(cyg_epsilon_y)*250)
        turtle.goto(float(cyg_zeta_x)*250, float(cyg_zeta_y)*250)
        turtle.goto(float(cyg_nu_x)*250, float(cyg_nu_y)*250)
        turtle.penup()
        turtle.goto(float(deneb_x)*250, float(deneb_y)*250)
        turtle.pendown()
        turtle.goto(float(cyg_xi_x)*250, float(cyg_xi_y)*250)


def hydra():
    with open('stars.txt', 'r') as original_star_file:
        hydra = original_star_file.readlines()
    hydra_stars = []

    for stars in hydra:
        if 'HYA ZETA' in stars:
            hydra_stars.append(stars)
        if 'HYA ETA' in stars:
            hydra_stars.append(stars)
        if 'HYA SIGMA' in stars:
            hydra_stars.append(stars)
        if 'HYA DELTA' in stars:
            hydra_stars.append(stars)
        if 'HYA EPSILON' in stars:
            hydra_stars.append(stars)
        if 'HYA THETA' in stars:
            hydra_stars.append(stars)
        if 'HYA IOTA' in stars:
            hydra_stars.append(stars)
        if 'HYA ALPHA' in stars:
            hydra_stars.append(stars)
        if 'HYA UPSIL1' in stars:
            hydra_stars.append(stars)
        if 'HYA UPSIL2' in stars:
            hydra_stars.append(stars)
        if 'HYA LAMBDA' in stars:
            hydra_stars.append(stars)
        if 'HYA MU' in stars:
            hydra_stars.append(stars)
        if 'HYA NU' in stars:
            hydra_stars.append(stars)
        if 'HYA BETA' in stars:
            hydra_stars.append(stars)
        if 'HYA XI' in stars:
            hydra_stars.append(stars)
        if 'HYA GAMMA' in stars:
            hydra_stars.append(stars)
        if 'HYA PI' in stars:
            hydra_stars.append(stars)

    with open('hydra_stars.txt', 'w') as hyd_stars:
        hyd_stars.writelines(hydra_stars)

    with open('hydra_stars.txt', 'r') as stars_hydra:
        for lines in stars_hydra:
            coordinates = lines.split()
            if 'HYA ZETA' in lines:
                hya_zeta_x = coordinates[0]
                hya_zeta_y = coordinates[1]
            if 'HYA ETA' in lines:
                hya_eta_x = coordinates[0]
                hya_eta_y = coordinates[1]
            if 'HYA SIGMA' in lines:
                hya_sigma_x = coordinates[0]
                hya_sigma_y = coordinates[1]
            if 'HYA DELTA' in lines:
                hya_delta_x = coordinates[0]
                hya_delta_y = coordinates[1]
            if 'HYA EPSILON' in lines:
                hya_epsilon_x = coordinates[0]
                hya_epsilon_y = coordinates[1]
            if 'HYA THETA' in lines:
                hya_theta_x = coordinates[0]
                hya_theta_y = coordinates[1]
            if 'HYA IOTA' in lines:
                hya_iota_x = coordinates[0]
                hya_iota_y = coordinates[1]
            if 'HYA ALPHA' in lines:
                hya_alpha_x = coordinates[0]
                hya_alpha_y = coordinates[1]
            if 'HYA UPSIL1' in lines:
                hya_upsil1_x = coordinates[0]
                hya_upsil1_y = coordinates[1]
            if 'HYA UPSIL2' in lines:
                hya_upsil2_x = coordinates[0]
                hya_upsil2_y = coordinates[1]
            if 'HYA LAMBDA' in lines:
                hya_lambda_x = coordinates[0]
                hya_lambda_y = coordinates[1]
            if 'HYA MU' in lines:
                hya_mu_x = coordinates[0]
                hya_mu_y = coordinates[1]
            if 'HYA NU' in lines:
                hya_nu_x = coordinates[0]
                hya_nu_y = coordinates[1]
            if 'HYA BETA' in lines:
                hya_beta_x = coordinates[0]
                hya_beta_y = coordinates[1]
            if 'HYA XI' in lines:
                hya_xi_x = coordinates[0]
                hya_xi_y = coordinates[1]
            if 'HYA GAMMA' in lines:
                hya_gamma_x = coordinates[0]
                hya_gamma_y = coordinates[1]
            if 'HYA PI' in lines:
                hya_pi_x = coordinates[0]
                hya_pi_y = coordinates[1]

        turtle.penup()
        turtle.goto(float(hya_zeta_x)*250, float(hya_zeta_y)*250)
        turtle.pendown()
        turtle.goto(float(hya_eta_x)*250, float(hya_eta_y)*250)
        turtle.goto(float(hya_sigma_x)*250, float(hya_sigma_y)*250)
        turtle.goto(float(hya_delta_x)*250, float(hya_delta_y)*250)
        turtle.goto(float(hya_epsilon_x)*250, float(hya_epsilon_y)*250)
        turtle.goto(float(hya_zeta_x)*250, float(hya_zeta_y)*250)
        turtle.goto(float(hya_theta_x)*250, float(hya_theta_y)*250)
        turtle.goto(float(hya_iota_x)*250, float(hya_iota_y)*250)
        turtle.goto(float(hya_alpha_x)*250, float(hya_alpha_y)*250)
        turtle.goto(float(hya_upsil1_x)*250, float(hya_upsil1_y)*250)
        turtle.goto(float(hya_upsil2_x)*250, float(hya_upsil2_y)*250)
        turtle.goto(float(hya_lambda_x)*250, float(hya_lambda_y)*250)
        turtle.goto(float(hya_mu_x)*250, float(hya_mu_y)*250)
        turtle.goto(float(hya_nu_x)*250, float(hya_nu_y)*250)
        turtle.goto(float(hya_alpha_x)*250, float(hya_alpha_y)*250)
        turtle.goto(float(hya_beta_x)*250, float(hya_beta_y)*250)
        turtle.goto(float(hya_xi_x)*250, float(hya_xi_y)*250)
        turtle.goto(float(hya_beta_x)*250, float(hya_beta_y)*250)
        turtle.goto(float(hya_gamma_x)*250, float(hya_gamma_y)*250)
        turtle.goto(float(hya_pi_x)*250, float(hya_pi_y)*250)


if __name__ == '__main__':
    main()
