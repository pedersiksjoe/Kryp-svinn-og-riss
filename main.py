from __future__ import print_function
import streamlit as st
import math
from mailmerge import MailMerge



#st.latex(r'''\tag*{hi} E_{cm}''')


tab1, tab2, tab3 = st.tabs(["Generell materialdata", "Kryp og svinn", "Rissberegning"])

tab3.header("Beregning av rissvidde")
#inndata
tab1.header("Geometri og laster")
h = tab1.number_input("Tversnittshøyde[mm]:",value=250)
b = tab1.number_input("Tversnittsbredde[mm]:",value=1000)
ø = tab1.number_input("Armeringsdiameter[mm]:",min_value=8, max_value=32, value=16)
s = tab1.number_input("Senteravstand armering[mm}:",value=200)
c = tab1.number_input("Nominell overdekning[mm]:",value=35)
M_ed_short = tab3.number_input("MEd_short[kNm]:",value=10.0)
M_ed_long = tab3.number_input("MEd_long[kNm]:",value=10.0)
M_ed_tot = M_ed_long+M_ed_short                     #Totalt moment
tab3.write("MEd_tot = "+str(M_ed_tot)+" kNm")

Ac = h*b                                            #Betongareal
As = (math.pow(ø/2,2)*math.pi*b)/s                  #Armeringsareal
d = h-c-(ø/2)                                       #Effektiv tverrsnittstykkelse
tab1.write("As = "+str(As))

RHstr = tab2.slider(
    "Relativ luftfuktighet (RH)",
    min_value=0,
    max_value=100,
    step=10,
    value=70
)
RH = int(RHstr)/100
RH0 = 1                             #(B.12)


t0 = tab2.number_input("Betongens alder ved pålasting (døgn):",value=28, step=1)

Levetid = tab2.radio(
    "Levetid",
    ("10", "25", "50", "100", "200", "300"),
    index=2
)
yr = 365
t = int(Levetid)*yr

sement = tab2.radio(
    "Sementklasse",
    ("S", "N", "R"),
    index=1
)

uttørking = tab2.radio(
    "Uttørking på 1 eller 2 sider",
    ("1 side", "2 sider"),
    index=1
)
if uttørking == "1 side":
    u = 1*b
elif uttørking == "2 sider":
    u = 2*b

tab3.header("Rissparametere iht 7.3.4(3):")
k_1 =tab3.number_input("k1 (0.8 for stenger med god heft):", value=0.8)
k_2 =tab3.number_input("k2 (0.5 for ren bøying, 1 for rent strekk):", value=0.5)


#Betongkvaliteter
Bx = ["f_ctm","E_cm","f_ck"]
B25 = [2.2,31,25]
B30 = [2.6,33,30]
B35 = [3.2,34,35]
B45 = [3.8,35,45]

tab1.header("Generelle materialdata:")
valg_betongkvalitet = tab1.radio(
    "Velg betongkvalitet",
    ("B25","B30","B35","B45")
)
if valg_betongkvalitet == "B25":
    betongkvalitet = B25
elif valg_betongkvalitet == "B30":
    betongkvalitet = B30
elif valg_betongkvalitet == "B35":
    betongkvalitet = B35
elif valg_betongkvalitet == "B45":
    betongkvalitet = B45

f_ctm = tab3.number_input("fctm[MPa]:",value=betongkvalitet[0])
E_cm = tab3.number_input("Ecm[GPa]:",value=betongkvalitet[1])
E_s = tab3.number_input("Es[GPa]:",value=200)
f_ck = betongkvalitet[2]
f_cm = f_ck+8



if uttørking == "1 side":                           #B.1(1)
    u = 1*b
elif uttørking == "2 sider":
    u = 2*b
h0 = (2*Ac)/(u)

#3.1.4 - Tabell 3.3
if h0 < 100:
    k_h = 1
elif h0 >= 100 and h0 < 200:
    k_h = (1-(h0-100)*0.0015)
elif h0 >= 200 and h0 < 300:
    k_h = (0.85-(h0-200)*0.001)
elif h0 >= 400 and h0 < 500:
    k_h = (0.75-(h0-300)*0.00025)
else:
    k_h = 0.7

#B.2
if sement == "S":
    alfa_ds1 = 3
    alfa_ds2 = 0.13
    alfa_t0 = -1
elif sement == "N":
    alfa_ds1 = 4
    alfa_ds2 = 0.12
    alfa_t0 = 0
elif sement == "R":
    alfa_ds1 = 6
    alfa_ds2 = 0.11
    alfa_t0 =1


f_cmo = 10                                          #(B.12)

beta_f_cm =16.8/(math.sqrt(f_cm))                   #(B.4)
beta_t0 = 1/(0.1+math.pow(t0,0.2))                  #(B.4)
alfa_1 = math.pow((35/f_cm),0.7)                    #(B.8.c)
alfa_2 = math.pow((35/f_cm),0.2)                    #(B.8.c)
alfa_3 = math.pow((35/f_cm),0.5)                    #(B.8.c)

#(B.3a)
#(B.3b)
# (B.8a)
# (B.8b)

if f_cm <= 35:
    phi_RH = 1+((1-RH)/(0.1*math.pow(h0,(1/3))))
    beta_H =min(1500,1.5*(1+math.pow(0.012*100*RH,18))*h0+250)
else:
    phi_RH = (1 + ((1 - RH) / (0.1 * math.pow(h0, (1 / 3))))*alfa_1)*alfa_2
    beta_H = min(1500*alfa_3, (1.5 * (1 + math.pow(0.012 * 100 * RH, 18)) * h0 + 250)*alfa_3)


beta_ct_to = math.pow((t-t0)/(beta_H+t-t0),0.3)                 #(B.7)
phi_0 = phi_RH*beta_f_cm*beta_t0                                #(B.2)
beta_RH = 1.55*(1-math.pow((RH/RH0),3))                         #(B.12)

#Nominell verdi for svinntøyning ved uttørking (B.11)
epsilon_cd0 = 0.85*((220+110*alfa_ds1)*math.pow(math.e,(-alfa_ds2*(f_cm/f_cmo))))*math.pow(10,-6)*beta_RH
beta_ds_t_ts = (t-t0)/((t-t0)+0.04*math.sqrt(math.pow(h0,3)))   #(3.10)
epsilon_cd_t = beta_ds_t_ts*k_h*epsilon_cd0                     #(3.9)
epsilon_ca_lim = 2.5*(f_ck-10)*math.pow(10,-6)                  #(3.12)
beta_ast = 1-math.pow(math.e,(-0.2*math.sqrt(t)))               #(3.15)
epsilon_ca_t = beta_ast*epsilon_ca_lim                          #(3.11)
epsilon_cs_beregnet = 1000*(epsilon_cd_t+epsilon_ca_t)          #(3.8)
phi_beregnet = phi_0*beta_ct_to                                 #(B.1)


phi = tab2.number_input("Kryptall:",value=round(phi_beregnet,3))
epsilon_cs = (tab2.number_input("Svinntøyning:",value=round(epsilon_cs_beregnet,3)))/1000000


n_short = M_ed_short/M_ed_tot                       #Andel korttidslast
n_long = M_ed_long/M_ed_tot                         #Andel langtidslast
E_c_short = E_cm                                    #E-modul korttidslast
E_c_long = E_cm/(1+phi)                             #E-modul langtidslast
eta_short = E_s/E_c_short
eta_long = E_s/E_c_long
E_c_mean = n_short*E_c_short + n_long*E_c_long      #E-modul mean
eta_mean = n_short*eta_short + n_long*eta_long


rho = As/(b*d)                                                  #Armeringstetthet
eta = eta_mean
rhoeta = rho*eta
alpha = math.sqrt(math.pow(rhoeta,2)+2*rhoeta)-rhoeta
I_s = As*(1-alpha)*(1-(alpha/3))*math.pow(d,2)
sigma_s = (M_ed_tot/I_s)*(1-alpha)*d                            #Spenning i armering fra moment



# Beregninger etter EC2
f_ct_eff = f_ctm                                                #Antatt >28 dagers herding
#7.3.4(2):
k_t = 0.6*n_long+0.4*n_short
#7.3.2(3):
h_c_ef1 = 2.5*(h-d)
h_c_ef2 = (h-alpha*d)/3
h_c_ef3 = h/2
h_c_ef = min(h_c_ef1, h_c_ef2, h_c_ef3)
A_c_eff = h_c_ef*b
#7.3.4(2):
alpha_e = E_s/E_cm
rho_p_eff = As/A_c_eff
#7.9:
epsilon_sm_e_cm1 = 0.6*(sigma_s/E_s)
epsilon_sm_e_cm2 = (sigma_s-k_t*(f_ct_eff/rho_p_eff)*(1+alpha_e*rho_p_eff))/E_s
epsilon_sm_e_cm = max(epsilon_sm_e_cm1, epsilon_sm_e_cm2)
#7.3.4(3):
k_3 = 3.4
k_4 = 0.425
#7.11
s_r_max1 = 1.3*(h-alpha*d)
s_r_max2 = k_3*c+(k_1*k_2*k_4*ø)/rho_p_eff
s_grenseverdi = 5*(c+(ø/2))

if s > s_grenseverdi:
    s_r_max = s_r_max1
else:
    s_r_max = s_r_max2

#(7.8)
w_k = 1000*s_r_max*epsilon_sm_e_cm

#med riss
w_k_ecs = 1000*s_r_max*(epsilon_sm_e_cm+epsilon_cs)

st.sidebar.header("Statiske beregninger")

st.sidebar.write("Armeringsspenning: "+str(round(sigma_s*1000000))+ " MPa")
st.sidebar.write("Rissvidde uten svinntøyning: "+str(round(w_k,3))+ " mm")
st.sidebar.write("Rissvidde med svinntøyning: "+str(round(w_k_ecs,3))+ " mm")

document = MailMerge("template.docx")
print(document.get_merge_fields())
dok = st.text_input("Dokumentnavn", value="Beregning av riss")
document.merge(
        h= str(h),
        b= str(b),
        A_c= str(Ac),
        d= str(round(d)),
        ø= str(ø),
        c_nom= str(c),
        Betongkvalitet= str(valg_betongkvalitet),
        E_c= str(E_cm),
        f_cm= str(f_cm),
        f_ctm= str(f_ctm),
        s= str(s),
        f_ck= str(f_ck)
    )
if st.button("Dokumentasjon"):
    file_temp = document.write("template_temp.docx")
    file = open("template_temp.docx", "rb")
    st.download_button(label="Last ned "+dok+".docx", file_name= dok+".docx", data=file, mime="application/vnd.openxmlformats-officedocument.wordprocessingml.document")
    file.close(