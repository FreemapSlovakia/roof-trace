# Budovy z leteckého laserového skenovania

Obrysy budov vektorizované z bodov triedy „budova“ leteckého laserového skenovania (LLS)
ÚGKK SR, pre doplnenie a spresnenie budov v OpenStreetMap. Obrysy idú po okraji strechy, nie po
múroch, a sú zarovnané: pravouhlé tam, kde je pravouhlá budova.

Susedné budovy sú rozdelené tam, kde sa ich strechy líšia (výškou, sklonom, smerom hrebeňa, úžľabím
medzi nimi). Delí sa skôr viac než menej: spojiť dve plochy v editore je rýchlejšie, než hľadať, kde
jednu rozdeliť.

Výstup vznikol programom [roof-trace](https://github.com/FreemapSlovakia/roof-trace) automaticky
a **nie je určený na priame nahratie do OSM**. Každú budovu treba pred nahratím skontrolovať
voči ortofotomape.

## Súbory

Územie je rozdelené na štvorce 20 × 20 km v S-JTSK (EPSG:8353). Adresár `<x>_<y>` je štvorec
s juhozápadným rohom v `x`, `y` a obsahuje:

- `new.osm` – budovy, ktoré v OSM ešte nie sú,
- `overlap.osm` – budovy, ktoré sa prekrývajú s existujúcou budovou v OSM, na porovnanie,
- vedľa adresára `<x>_<y>.geojson` – všetky obrysy štvorca v S-JTSK, aj s diagnostickými
  atribútmi.

`done.geojson` (WGS 84) ukazuje hotové štvorce; atribút `dir` je názov adresára. Ďalšie pribúdajú,
kým beží spracovanie celého Slovenska.

Budova patrí do štvorca, v ktorom leží jej ťažisko, takže sa susedné štvorce neprekrývajú. Budovy
na rozhraní štvorcov nezdieľajú body so susedným štvorcom.

## Atribúty

- `building=yes`,
- `height` – výška hrebeňa nad najnižším terénom okolo budovy (m),
- `roof:shape`, `roof:height`, `roof:orientation` – kde je strecha jedným z jednoduchých tvarov,
- `source=LLS: ÚGKK SR` – geometria aj všetky atribúty sú z LLS.

Rohy susedných budov bližšie ako 0,4 m sú jeden bod a roh na stene susednej budovy je bodom aj tej
steny, takže pristavané budovy majú spoločné steny.

## Použitie v JOSM

Súbory majú `upload="never"`, JOSM ich teda nedovolí nahrať priamo. Vybrané budovy skopírujte
(Ctrl+C) do vrstvy s dátami z OSM a vložte (Ctrl+Shift+V), skontrolujte a až potom nahrajte.

Údaje LLS sú z rôznych rokov. Budova môže byť medzičasom zbúraná alebo prestavaná a novšie budovy
chýbajú.

## Licencia

Zdroj údajov: ÚGKK SR, letecké laserové skenovanie. Pri použití uveďte „LLS: ÚGKK SR“.
